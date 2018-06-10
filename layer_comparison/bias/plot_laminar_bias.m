function plot_laminar_bias(subjects, contrast, varargin)

% Parse inputs
defaults = struct('data_dir','d:/meg_laminar/derivatives/spm12',...
    'surf_dir', 'D:/meg_laminar/derivatives/freesurfer','inv_type','EBB',...
    'patch_size',0.4,'thresh_percentile',80,'roi_type','mean',...
    'whole_brain', false, 'plot_ext','',...
    'recompute_roi',false, 'ylim', [], 'perc', false, 'perc_limit', 0);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');
addpath('D:\meg_laminar\layer_comparison');

thresh_type='lower';
switch contrast.comparison_name
    case 'dots_beta_erd'
        thresh_type='upper';
    case 'dots_alpha'
        thresh_type='upper';
end

depth_pial_wm_diff_cv=[];
lfn_pial_wm_diff_cv=[];
all_lfn_diff=[];
all_depth_lfn_diff=[];
all_depth_diff=[];
all_lfn_depth_diff=[];
all_depth_subj=[];
all_lfn_subj=[];

for s_idx=1:length(subjects)
    subj_info=subjects(s_idx);
    subj_surf_dir=fullfile(params.surf_dir,subj_info.subj_id);
    
    orig_white_mesh=fullfile(subj_surf_dir,'white.hires.deformed.surf.gii');
    white_mesh=fullfile(subj_surf_dir,'ds_white.hires.deformed.surf.gii');
    white_inflated=fullfile(subj_surf_dir,'ds_white.hires.deformed_inflated.surf.gii');
    orig_pial_mesh=fullfile(subj_surf_dir,'pial.hires.deformed.surf.gii');
    pial_mesh=fullfile(subj_surf_dir,'ds_pial.hires.deformed.surf.gii');
    pial_inflated=fullfile(subj_surf_dir,'ds_pial.hires.deformed_inflated.surf.gii');
    pial_white_map=map_pial_to_white(white_mesh, pial_mesh, 'mapType', 'link',...
        'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh);
    white_pial_map=map_white_to_pial(white_mesh, pial_mesh, 'mapType', 'link',...
        'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh);
    pial_hemisphere_map=get_hemisphere_map(pial_mesh, orig_pial_mesh);
    white_hemisphere_map=get_hemisphere_map(white_mesh, orig_white_mesh);
    pial=gifti(pial_mesh);
    wm=gifti(white_mesh);
    nvertices=size(pial.vertices,1);

    first_session=1;
    if strcmp(subj_info.subj_id,'nc')
        first_session=3;
    end
    foi_dir=fullfile(params.data_dir, subj_info.subj_id,...
        sprintf('ses-%02d',first_session), 'grey_coreg', params.inv_type,....
        ['p' num2str(params.patch_size)], contrast.zero_event,...
        ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);

    lfn_filename=fullfile(foi_dir, sprintf('br%s_%d.mat',subj_info.subj_id, first_session));    
    grey_lfn_fname=fullfile(subj_surf_dir, 'ds_grey_lfn.mat');
    
    if exist(grey_lfn_fname,'file')~=2
        D=spm_eeg_load(lfn_filename);
        lfn=sqrt(sum(D.inv{1}.inverse.L.^2,1))';        
        save(grey_lfn_fname,'lfn');
    else
        load(grey_lfn_fname);
    end
    pial_lfn=lfn(nvertices+1:end);
    pial_lfn=spm_mesh_smooth(pial, pial_lfn, 20);
    wm_lfn=lfn(1:nvertices);
    wm_lfn=spm_mesh_smooth(wm, wm_lfn, 20);
    mapped_wm_lfn=wm_lfn(pial_white_map);
    lfn_diff=pial_lfn-mapped_wm_lfn;
    if params.perc
        lfn_diff=lfn_diff./mean([pial_lfn mapped_wm_lfn],2).*100.0;
    end

    pial_metric_mask=find(pial_lfn>22);
    mapped_wm_metric_mask=find(mapped_wm_lfn>19);
    metric_mask=union(pial_metric_mask,mapped_wm_metric_mask);
    pial_white_dists=sqrt(sum((pial.vertices-wm.vertices(pial_white_map,:)).^2,2));
    metric_mask=intersect(metric_mask, find(pial_white_dists>=0.11));


    hemisphere=contrast.hemisphere;
    region=contrast.region;
    if params.whole_brain
        hemisphere='';
        region='';
    end
    foi_dir=fullfile(params.data_dir, subj_info.subj_id,...
            'grey_coreg', params.inv_type,....
            ['p' num2str(params.patch_size)], contrast.zero_event,...
            ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
    [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, contrast.comparison_name, ...
            thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
            pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',params.thresh_percentile,...
            'type',params.roi_type, 'region', region, 'hemisphere', hemisphere,...
            'pial_hemisphere_map', pial_hemisphere_map,...
            'white_hemisphere_map', white_hemisphere_map, 'recompute', params.recompute_roi);
    mask=intersect(metric_mask,mask);
       
    pial_wm_diff=gifti(fullfile(foi_dir,['pial-white.' contrast.comparison_name '.diff.gii']));
    pial_wm_diff_data=pial_wm_diff.cdata(:,:);
    pial_wm_diff_cv=mean(pial_wm_diff_data,2)./std(pial_wm_diff_data,[],2);

    pial_depth_fname=fullfile(subj_surf_dir,'ds_pial_sulcal_depth.mat');
    wm_depth_fname=fullfile(subj_surf_dir,'ds_white_sulcal_depth.mat');
    if exist(pial_depth_fname,'file')~=1 || exist(wm_depth_fname,'file')~=2
        [pial_depth,HS]=compute_sulcal_depth(pial_mesh);
        mapping=dsearchn(HS.vertices,wm.vertices);
        wm_depth=sqrt(sum((wm.vertices-HS.vertices(mapping,:)).^2,2));
        save(pial_depth_fname,'pial_depth');
        save(wm_depth_fname,'wm_depth');
    else
        load(pial_depth_fname);
        load(wm_depth_fname);
    end
    mapped_wm_depth=wm_depth(pial_white_map);
    depth_diff=pial_depth-mapped_wm_depth;
    if params.perc
        depth_diff=depth_diff./mean([pial_depth mapped_wm_depth],2).*100.0;
    end
    
    % Add depth and LFN to mask
    if params.perc && params.perc_limit>0
        depth_mask=intersect(mask, find(abs(depth_diff)<=params.perc_limit));    
        depth_pial_wm_diff_cv(end+1:end+length(depth_mask))=pial_wm_diff_cv(depth_mask);
        all_depth_diff(end+1:end+length(depth_mask))=depth_diff(depth_mask);
        all_depth_lfn_diff(end+1:end+length(depth_mask))=lfn_diff(depth_mask);
        all_depth_subj(end+1:end+length(depth_mask))=s_idx.*ones(length(depth_mask),1);
        
        lfn_mask=intersect(mask, find(abs(lfn_diff)<=params.perc_limit));    
        lfn_pial_wm_diff_cv(end+1:end+length(lfn_mask))=pial_wm_diff_cv(lfn_mask);
        all_lfn_diff(end+1:end+length(lfn_mask))=lfn_diff(lfn_mask);
        all_lfn_depth_diff(end+1:end+length(lfn_mask))=depth_diff(lfn_mask);
        all_lfn_subj(end+1:end+length(lfn_mask))=s_idx.*ones(length(lfn_mask),1);      
    else
        depth_pial_wm_diff_cv(end+1:end+length(mask))=pial_wm_diff_cv(mask);
        all_depth_diff(end+1:end+length(mask))=depth_diff(mask);
        all_depth_lfn_diff(end+1:end+length(mask))=lfn_diff(mask);
        all_depth_subj(end+1:end+length(mask))=s_idx.*ones(length(mask),1);
        
        lfn_pial_wm_diff_cv(end+1:end+length(mask))=pial_wm_diff_cv(mask);
        all_lfn_diff(end+1:end+length(mask))=lfn_diff(mask);
        all_lfn_depth_diff(end+1:end+length(mask))=depth_diff(mask);
        all_lfn_subj(end+1:end+length(mask))=s_idx.*ones(length(mask),1);      
    end       
    
end

colors=[102 102 102;27 158 119;217 95 2;117 112 179;231 41 138;102 166 30;230 171 2;166 118 29]./255;


%%
% LFN diff - laminar pref scatterhist
fig=figure();
% h=scatterhist(all_lfn_diff,lfn_pial_wm_diff_cv,'Group',all_lfn_subj,...
%     'Kernel','on','Location','SouthEast','Direction','out','Marker','.',...
%     'MarkerSize',3,'Color',colors,'LineStyle',{'-'},'LineWidth',1);
xlims=[];
ylims=[];
if params.perc && params.perc_limit>0
    xlims=[-params.perc_limit params.perc_limit];
end
if length(params.ylim)>0
    ylims=params.ylim;    
end
h=myscatterhist(all_lfn_diff,lfn_pial_wm_diff_cv,'Group',all_lfn_subj,...
    'Kernel','on','Location','SouthEast','Direction','out','Marker','.',...
    'MarkerSize',3,'Color',colors,'LineStyle',{'-'},'LineWidth',1,'Xlim',xlims,'Ylim',ylims);
axes(h(1));
hold all;
subj_r=[];
subj_z=[];
subj_partial_r=[];
subj_partial_z=[];
for s_idx=1:length(subjects)
    idx=find(all_lfn_subj==s_idx);
    
    [rho,p]=corr(all_lfn_diff(idx)',lfn_pial_wm_diff_cv(idx)','Type','Spearman');
    subj_r(s_idx)=rho;
    subj_z(s_idx)=0.5*log((1+rho)/(1-rho));
    disp(sprintf('LFN diff - Pial/WM diff, r=%.4f, p=%.4f\n', rho, p));
    
    [rho,p]=partialcorr(all_lfn_diff(idx)',lfn_pial_wm_diff_cv(idx)',all_lfn_depth_diff(idx)','Type','Spearman');
    subj_partial_r(s_idx)=rho;
    subj_partial_z(s_idx)=0.5*log((1+rho)/(1-rho));
    disp(sprintf('LFN diff - Pial/WM diff, partial r=%.4f, p=%.4f\n', rho, p));
    
    pPoly = polyfit(all_lfn_diff(idx),lfn_pial_wm_diff_cv(idx), 1); % Linear fit of xdata vs ydata
    linePointsX = [min(all_lfn_diff(idx)) max(all_lfn_diff(idx))]; % find left and right x values
    linePointsY = polyval(pPoly,linePointsX); 
    plot (linePointsX,linePointsY,'--','Color',colors(s_idx,:),'LineWidth',2);
end
main_xl=xlim();
main_yl=ylim();
[hyp,p,ci,stats]=ttest(subj_z);
disp(sprintf('All LFN diff - Pial/WM diff, mean r=%.4f, t(%d)=%.4f, p=%.4f\n', mean(subj_r), stats.df, stats.tstat, p));
[hyp,p,ci,stats]=ttest(subj_partial_z);
disp(sprintf('All LFN diff - Pial/WM diff, mean partial r=%.4f, t(%d)=%.4f, p=%.4f\n', mean(subj_partial_r), stats.df, stats.tstat, p));
axes(h(2));
hold all;
yl=ylim();
for s_idx=1:length(subjects)
    idx=find(all_lfn_subj==s_idx);
    subj_mean=median(all_lfn_diff(idx));
    plot([subj_mean subj_mean],yl,'--','Color',colors(s_idx,:),'LineWidth',2);
end
xlim(main_xl);
axes(h(3));
hold all;
yl=ylim();
for s_idx=1:length(subjects)
    idx=find(all_lfn_subj==s_idx);
    subj_mean=median(lfn_pial_wm_diff_cv(idx));
    plot([subj_mean subj_mean],yl,'--','Color',colors(s_idx,:),'LineWidth',2);
end
xlim(main_yl);
xlabel('Pial LFN - WM LFN','FontSize',14,'FontName','Arial');
ylabel('Effect size(|Pial Power D| - |WM Power D|)','FontSize',14,'FontName','Arial');
set(gca,'FontName','Arial');

%%
% Depth diff - laminar pref scatterhist
fig=figure();
xlims=[];
ylims=[];
if params.perc && params.perc_limit>0
    xlims=[-params.perc_limit params.perc_limit];
end
if length(params.ylim)>0
    ylims=params.ylim;    
end
h=myscatterhist(all_depth_diff,depth_pial_wm_diff_cv,'Group',all_depth_subj,...
    'Kernel','on','Location','SouthEast','Direction','out','Marker','.',...
    'MarkerSize',3,'Color',colors,'LineStyle',{'-'},'LineWidth',1,'Xlim',xlims,'Ylim',ylims);
axes(h(1));
hold all;
subj_r=[];
subj_z=[];
subj_partial_r=[];
subj_partial_z=[];
for s_idx=1:length(subjects)
    idx=find(all_depth_subj==s_idx);
    
    [rho,p]=corr(all_depth_diff(idx)',depth_pial_wm_diff_cv(idx)','Type','Spearman');
    subj_r(s_idx)=rho;
    subj_z(s_idx)=0.5*log((1+rho)/(1-rho));
    disp(sprintf('Depth diff - Pial/WM diff, r=%.4f, p=%.4f\n', rho, p));
    
    [rho,p]=partialcorr(all_depth_diff(idx)',depth_pial_wm_diff_cv(idx)',all_depth_lfn_diff(idx)','Type','Spearman');
    subj_partial_r(s_idx)=rho;
    subj_partial_z(s_idx)=0.5*log((1+rho)/(1-rho));
    disp(sprintf('Depth diff - Pial/WM diff, partial r=%.4f, p=%.4f\n', rho, p));
    
    pPoly = polyfit(all_depth_diff(idx),depth_pial_wm_diff_cv(idx), 1); % Linear fit of xdata vs ydata
    linePointsX = [min(all_depth_diff(idx)) max(all_depth_diff(idx))]; % find left and right x values
    linePointsY = polyval(pPoly,linePointsX); 
    plot(linePointsX,linePointsY,'--','Color',colors(s_idx,:),'LineWidth',2);
end
main_xl=xlim();
main_yl=ylim();
[hyp,p,ci,stats]=ttest(subj_z);
disp(sprintf('All Depth diff - Pial/WM diff, mean r=%.4f, t(%d)=%.4f, p=%.4f\n', mean(subj_r), stats.df, stats.tstat, p));
[hyp,p,ci,stats]=ttest(subj_partial_z);
disp(sprintf('All Depth diff - Pial/WM diff, mean partial r=%.4f, t(%d)=%.4f, p=%.4f\n', mean(subj_partial_r), stats.df, stats.tstat, p));
axes(h(2));
hold all;
yl=ylim();
for s_idx=1:length(subjects)
    idx=find(all_depth_subj==s_idx);
    subj_mean=median(all_depth_diff(idx));
    plot([subj_mean subj_mean],yl,'--','Color',colors(s_idx,:),'LineWidth',2);
end
xlim(main_xl);
axes(h(3));
hold all;
yl=ylim();
for s_idx=1:length(subjects)
    idx=find(all_depth_subj==s_idx);
    subj_mean=median(depth_pial_wm_diff_cv(idx));
    plot([subj_mean subj_mean],yl,'--','Color',colors(s_idx,:),'LineWidth',2);
end
xlim(main_yl);
xlabel('Pial depth - WM depth','FontSize',14,'FontName','Arial');
ylabel('CSV(|Pial Power D| - |WM Power D|)','FontSize',14,'FontName','Arial');
set(gca,'FontName','Arial');

rmpath('D:\meg_laminar\layer_comparison');