function plot_laminar_bias(subjects, contrast, sub_dir, varargin)

% Parse inputs
defaults = struct('data_dir','d:/pred_coding',...
    'surf_dir', 'D:/pred_coding/surf','inv_type','EBB',...
    'patch_size',0.4,'thresh_percentile',80,'roi_type','mean',...
    'whole_brain', false, 'plot_ext','',...
    'recompute_roi',false, 'ylim', []);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');
addpath('D:\pred_coding\src\matlab\analysis\layer_comparison');

thresh_type='lower';
switch contrast.comparison_name
    case 'dots_beta_erd'
        thresh_type='upper';
    case 'dots_alpha'
        thresh_type='upper';
end

plot_dir=fullfile('C:\Users\jbonai\Dropbox\meg\pred_coding\plots\layer_comparison',...
    contrast.comparison_name,'bias',sub_dir);
mkdir(plot_dir);

all_pial_wm_diff_cv=[];
all_lfn_diff=[];
all_depth_diff=[];
all_subj=[];

fid=fopen(fullfile(plot_dir,'stats.txt'),'w');

for s_idx=1:length(subjects)
    subj_info=subjects(s_idx);
    subj_surf_dir=fullfile(params.surf_dir,sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf');
    
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
    foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
        num2str(first_session), 'grey_coreg', params.inv_type,....
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

    pial_metric_mask=find(pial_lfn>22);
    mapped_wm_metric_mask=find(mapped_wm_lfn>19);
    metric_mask=union(pial_metric_mask,mapped_wm_metric_mask);

    hemisphere=contrast.hemisphere;
    region=contrast.region;
    if params.whole_brain
        hemisphere='';
        region='';
    end
    foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
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
    
    fprintf(fid,sprintf('Subject %d: %s - mean pial LFN -white LFN=%.4f\n', s_idx, contrast.comparison_name, full(mean(lfn_diff(mask)))));
    fprintf(fid,sprintf('Subject %d: %s - mean pial depth -white depth=%.4f\n', s_idx, contrast.comparison_name, mean(depth_diff(mask))));
   
    all_pial_wm_diff_cv(end+1:end+length(mask))=pial_wm_diff_cv(mask);
    all_lfn_diff(end+1:end+length(mask))=lfn_diff(mask);
    all_depth_diff(end+1:end+length(mask))=depth_diff(mask);
    all_subj(end+1:end+length(mask))=s_idx.*ones(length(mask),1);
end

colors=[102 102 102;27 158 119;217 95 2;117 112 179;231 41 138;102 166 30;230 171 2;166 118 29]./255;


%%
% Depth diff - LFN diff scatterhist
%
fig=figure();
h=scatterhist(all_depth_diff,all_lfn_diff,'Group',all_subj,...
    'Kernel','on','Location','SouthEast','Direction','out','Marker','.',...
    'MarkerSize',3,'Color',colors,'LineStyle',{'-'},'LineWidth',1);
axes(h(1));
main_xl=xlim();
main_yl=ylim();
hold all;
subj_r=[];
subj_z=[];
for s_idx=1:length(subjects)
    idx=find(all_subj==s_idx);
    [rho,p]=corr(all_depth_diff(idx)',all_lfn_diff(idx)','Type','Spearman');
    subj_r(s_idx)=rho;
    subj_z(s_idx)=0.5*log((1+rho)/(1-rho));
    fprintf(fid,sprintf('Depth diff - LFN diff, r=%.4f, p=%.4f\n', rho, p));
    pPoly = polyfit(all_depth_diff(idx), all_lfn_diff(idx),1); % Linear fit of xdata vs ydata
    linePointsX = [min(all_depth_diff(idx)) max(all_depth_diff(idx))]; % find left and right x values
    linePointsY = polyval(pPoly,linePointsX); 
    plot (linePointsX,linePointsY,'--','Color',colors(s_idx,:),'LineWidth',2);
end
[hyp,p,ci,stats]=ttest(subj_z);
fprintf(fid,sprintf('All depth diff - LFN diff, mean r=%.4f, t(%d)=%.4f, p=%.4f\n', mean(subj_r), stats.df, stats.tstat, p));  
axes(h(2));
hold all;
yl=ylim();
for s_idx=1:length(subjects)
    idx=find(all_subj==s_idx);
    subj_mean=median(all_depth_diff(idx));
    plot([subj_mean subj_mean],yl,'--','Color',colors(s_idx,:),'LineWidth',2);
end
axes(h(3));
hold all;
yl=ylim();
for s_idx=1:length(subjects)
    idx=find(all_subj==s_idx);
    subj_mean=median(all_lfn_diff(idx));
    plot([subj_mean subj_mean],yl,'--','Color',colors(s_idx,:),'LineWidth',2);
end
xlabel('Pial depth - WM depth','FontSize',14,'FontName','Arial');
ylabel('Pial LFN - WM LFN','FontSize',14,'FontName','Arial');
set(gca,'FontName','Arial');
saveas(fig, fullfile(plot_dir, 'lfn_diff_depth_diff.eps'), 'epsc');

%%
% Depth diff - LFN diff bivariate histogram
fig=figure();
nbins=100;
depth_bin_size=(main_xl(2)-main_xl(1))/nbins;
lfn_bin_size=(main_yl(2)-main_yl(1))/nbins;
X=[all_depth_diff', all_lfn_diff'];
depth_ctrs=[main_xl(1)+depth_bin_size/2:depth_bin_size:main_xl(2)-depth_bin_size/2];
lfn_ctrs=[main_yl(1)+lfn_bin_size/2:lfn_bin_size:main_yl(2)-lfn_bin_size/2];
allN=[];
for s_idx=1:length(subjects)
    idx=find(all_subj==s_idx);
    [N,c]=hist3(X(idx,:),'Ctrs',{depth_ctrs lfn_ctrs});
    x_diff=diff(depth_ctrs);
    y_diff=diff(lfn_ctrs);
    x = repmat([x_diff, x_diff(end)], length(lfn_ctrs),1)';
    y = repmat([y_diff, y_diff(end)], length(depth_ctrs),1);
    % volume of the histogram
    V_tot  = sum(sum(x.*y.*N));
    allN(s_idx,:,:)=N./V_tot.*100.0;
end
imagesc(depth_ctrs,lfn_ctrs,squeeze(mean(allN))');
set(gca,'ydir','normal');
xlim(main_xl);
ylim(main_yl);
colorbar('westoutside');
xlabel('Pial depth - WM depth','FontSize',14,'FontName','Arial');
ylabel('Pial LFN - WM LFN','FontSize',14,'FontName','Arial');
set(gca,'FontName','Arial');
saveas(fig, fullfile(plot_dir, 'lfn_diff_depth_diff_hist.eps'), 'epsc');
saveas(fig, fullfile(plot_dir, 'lfn_diff_depth_diff_hist.png'), 'png');

   
%%
% LFN diff - laminar pref scatterhist
fig=figure();
h=scatterhist(all_lfn_diff,all_pial_wm_diff_cv,'Group',all_subj,...
    'Kernel','on','Location','SouthEast','Direction','out','Marker','.',...
    'MarkerSize',3,'Color',colors,'LineStyle',{'-'},'LineWidth',1);
axes(h(1));
hold all;
subj_r=[];
subj_z=[];
subj_partial_r=[];
subj_partial_z=[];
for s_idx=1:length(subjects)
    idx=find(all_subj==s_idx);
    
    [rho,p]=corr(all_lfn_diff(idx)',all_pial_wm_diff_cv(idx)','Type','Spearman');
    subj_r(s_idx)=rho;
    subj_z(s_idx)=0.5*log((1+rho)/(1-rho));
    fprintf(fid,sprintf('LFN diff - Pial/WM diff, r=%.4f, p=%.4f\n', rho, p));
    
    [rho,p]=partialcorr(all_lfn_diff(idx)',all_pial_wm_diff_cv(idx)',all_depth_diff(idx)','Type','Spearman');
    subj_partial_r(s_idx)=rho;
    subj_partial_z(s_idx)=0.5*log((1+rho)/(1-rho));
    fprintf(fid,sprintf('LFN diff - Pial/WM diff, partial r=%.4f, p=%.4f\n', rho, p));
    
    pPoly = polyfit(all_lfn_diff(idx),all_pial_wm_diff_cv(idx), 1); % Linear fit of xdata vs ydata
    linePointsX = [min(all_lfn_diff(idx)) max(all_lfn_diff(idx))]; % find left and right x values
    linePointsY = polyval(pPoly,linePointsX); 
    plot (linePointsX,linePointsY,'--','Color',colors(s_idx,:),'LineWidth',2);
end
if length(params.ylim)>0
    ylim(params.ylim);    
end
main_xl=xlim();
main_yl=ylim();
[hyp,p,ci,stats]=ttest(subj_z);
fprintf(fid,sprintf('All LFN diff - Pial/WM diff, mean r=%.4f, t(%d)=%.4f, p=%.4f\n', mean(subj_r), stats.df, stats.tstat, p));
[hyp,p,ci,stats]=ttest(subj_partial_z);
fprintf(fid,sprintf('All LFN diff - Pial/WM diff, mean partial r=%.4f, t(%d)=%.4f, p=%.4f\n', mean(subj_partial_r), stats.df, stats.tstat, p));
axes(h(2));
hold all;
yl=ylim();
for s_idx=1:length(subjects)
    idx=find(all_subj==s_idx);
    subj_mean=median(all_lfn_diff(idx));
    plot([subj_mean subj_mean],yl,'--','Color',colors(s_idx,:),'LineWidth',2);
end
axes(h(3));
hold all;
yl=ylim();
for s_idx=1:length(subjects)
    idx=find(all_subj==s_idx);
    subj_mean=median(all_pial_wm_diff_cv(idx));
    plot([subj_mean subj_mean],yl,'--','Color',colors(s_idx,:),'LineWidth',2);
end
if length(params.ylim)>0
    xlim(params.ylim);
end
xlabel('Pial LFN - WM LFN','FontSize',14,'FontName','Arial');
ylabel('CV(|Pial Power D| - |WM Power D|)','FontSize',14,'FontName','Arial');
set(gca,'FontName','Arial');
saveas(fig, fullfile(plot_dir, 'lfn_diff.eps'), 'epsc');

%%
% LFN diff - laminar pref bivariate histogram
fig=figure();
nbins=100;
lfn_bin_size=(main_xl(2)-main_xl(1))/nbins;
pial_wm_bin_size=(main_yl(2)-main_yl(1))/nbins;
X=[all_lfn_diff', all_pial_wm_diff_cv'];
lfn_ctrs=[main_xl(1)+lfn_bin_size/2:lfn_bin_size:main_xl(2)-lfn_bin_size/2];
pial_wm_ctrs=[main_yl(1)+pial_wm_bin_size/2:pial_wm_bin_size:main_yl(2)-pial_wm_bin_size/2];
allN=[];
for s_idx=1:length(subjects)
    idx=find(all_subj==s_idx);
    [N,c]=hist3(X(idx,:),'Ctrs',{lfn_ctrs pial_wm_ctrs});
    x_diff=diff(lfn_ctrs);
    y_diff=diff(pial_wm_ctrs);
    x = repmat([x_diff, x_diff(end)], length(lfn_ctrs),1)';
    y = repmat([y_diff, y_diff(end)], length(pial_wm_ctrs),1);
    % volume of the histogram
    V_tot  = sum(sum(x.*y.*N));
    allN(s_idx,:,:)=N./V_tot.*100.0;
end
imagesc(lfn_ctrs,pial_wm_ctrs,squeeze(mean(allN))');
set(gca,'ydir','normal');
xlim(main_xl);
ylim(main_yl);
colorbar('westoutside');
xlabel('Pial LFN - WM LFN','FontSize',14,'FontName','Arial');
ylabel('CV(|Pial Power D| - |WM Power D|)','FontSize',14,'FontName','Arial');
set(gca,'FontName','Arial');
saveas(fig, fullfile(plot_dir, 'lfn_diff_hist.eps'), 'epsc');
saveas(fig, fullfile(plot_dir, 'lfn_diff_hist.png'), 'png');

%%
% Depth diff - laminar pref scatterhist
fig=figure();
h=scatterhist(all_depth_diff,all_pial_wm_diff_cv,'Group',all_subj,...
    'Kernel','on','Location','SouthEast','Direction','out','Marker','.',...
    'MarkerSize',3,'Color',colors,'LineStyle',{'-'},'LineWidth',1);
axes(h(1));
hold all;
subj_r=[];
subj_z=[];
subj_partial_r=[];
subj_partial_z=[];
for s_idx=1:length(subjects)
    idx=find(all_subj==s_idx);
    
    [rho,p]=corr(all_depth_diff(idx)',all_pial_wm_diff_cv(idx)','Type','Spearman');
    subj_r(s_idx)=rho;
    subj_z(s_idx)=0.5*log((1+rho)/(1-rho));
    fprintf(fid,sprintf('Depth diff - Pial/WM diff, r=%.4f, p=%.4f\n', rho, p));
    
    [rho,p]=partialcorr(all_depth_diff(idx)',all_pial_wm_diff_cv(idx)',all_lfn_diff(idx)','Type','Spearman');
    subj_partial_r(s_idx)=rho;
    subj_partial_z(s_idx)=0.5*log((1+rho)/(1-rho));
    fprintf(fid,sprintf('Depth diff - Pial/WM diff, partial r=%.4f, p=%.4f\n', rho, p));
    
    pPoly = polyfit(all_depth_diff(idx),all_pial_wm_diff_cv(idx), 1); % Linear fit of xdata vs ydata
    linePointsX = [min(all_depth_diff(idx)) max(all_depth_diff(idx))]; % find left and right x values
    linePointsY = polyval(pPoly,linePointsX); 
    plot(linePointsX,linePointsY,'--','Color',colors(s_idx,:),'LineWidth',2);
end
if length(params.ylim)>0
    ylim(params.ylim);
end
main_xl=xlim();
main_yl=ylim();
[hyp,p,ci,stats]=ttest(subj_z);
fprintf(fid,sprintf('All Depth diff - Pial/WM diff, mean r=%.4f, t(%d)=%.4f, p=%.4f\n', mean(subj_r), stats.df, stats.tstat, p));
[hyp,p,ci,stats]=ttest(subj_partial_z);
fprintf(fid,sprintf('All Depth diff - Pial/WM diff, mean partial r=%.4f, t(%d)=%.4f, p=%.4f\n', mean(subj_partial_r), stats.df, stats.tstat, p));
axes(h(2));
hold all;
yl=ylim();
for s_idx=1:length(subjects)
    idx=find(all_subj==s_idx);
    subj_mean=median(all_depth_diff(idx));
    plot([subj_mean subj_mean],yl,'--','Color',colors(s_idx,:),'LineWidth',2);
end
axes(h(3));
hold all;
yl=ylim();
for s_idx=1:length(subjects)
    idx=find(all_subj==s_idx);
    subj_mean=median(all_pial_wm_diff_cv(idx));
    plot([subj_mean subj_mean],yl,'--','Color',colors(s_idx,:),'LineWidth',2);
end
if length(params.ylim)>0
    xlim(params.ylim);
end
xlabel('Pial depth - WM depth','FontSize',14,'FontName','Arial');
ylabel('CSV(|Pial Power D| - |WM Power D|)','FontSize',14,'FontName','Arial');
set(gca,'FontName','Arial');
saveas(fig, fullfile(plot_dir, 'depth_diff.eps'), 'epsc');

%%
% Depth diff - laminar pref bivariate histogram
fig=figure();
nbins=100;
depth_bin_size=(main_xl(2)-main_xl(1))/nbins;
pial_wm_bin_size=(main_yl(2)-main_yl(1))/nbins;
X=[all_depth_diff', all_pial_wm_diff_cv'];
depth_ctrs=[main_xl(1)+depth_bin_size/2:depth_bin_size:main_xl(2)-depth_bin_size/2];
pial_wm_ctrs=[main_yl(1)+pial_wm_bin_size/2:pial_wm_bin_size:main_yl(2)-pial_wm_bin_size/2];
allN=[];
for s_idx=1:length(subjects)
    idx=find(all_subj==s_idx);
    [N,c]=hist3(X(idx,:),'Ctrs',{depth_ctrs pial_wm_ctrs});
    x_diff=diff(depth_ctrs);
    y_diff=diff(pial_wm_ctrs);
    x = repmat([x_diff, x_diff(end)], length(lfn_ctrs),1)';
    y = repmat([y_diff, y_diff(end)], length(pial_wm_ctrs),1);
    % volume of the histogram
    V_tot  = sum(sum(x.*y.*N));
    allN(s_idx,:,:)=N./V_tot.*100.0;
end
imagesc(depth_ctrs,pial_wm_ctrs,squeeze(mean(allN))');
set(gca,'ydir','normal');
xlim(main_xl);
ylim(main_yl);
colorbar('westoutside');
xlabel('Pial depth - WM depth','FontSize',14,'FontName','Arial');
ylabel('CV(|Pial Power D| - |WM Power D|)','FontSize',14,'FontName','Arial');
set(gca,'FontName','Arial');
saveas(fig, fullfile(plot_dir, 'depth_diff_hist.eps'), 'epsc');
saveas(fig, fullfile(plot_dir, 'depth_diff_hist.png'), 'png');

rmpath('D:\pred_coding\src\matlab\analysis\layer_comparison');