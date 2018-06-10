function [depth_intercepts,lfn_intercepts,depth_means,lfn_means]=plot_subject_laminar_bias(subj_info, contrasts, varargin)

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

contrast_order=[3 4 6 5 1 2];
colors=[102,194,165;252,141,98;141,160,203;231,138,195;166,216,84;255,217,47]./255;
styles={'-','-','-','--','--','--'};

all_lfn_pial_wm_diff_data=[];
all_depth_pial_wm_diff_data=[];
all_lfn_diff=[];
all_depth_diff=[];
all_lfn_contrast=[];
all_depth_contrast=[];

subj_surf_dir=fullfile(params.surf_dir, subj_info.subj_id);

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
        
for c_idx=1:length(contrast_order)
    contrast=contrasts(contrast_order(c_idx));

    thresh_type='lower';
    switch contrast.comparison_name
        case 'dots_beta_erd'
            thresh_type='upper';
        case 'dots_alpha'
            thresh_type='upper';
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
        lfn_diff=lfn_diff./pial_lfn.*100.0;
    end

    pial_metric_mask=find(pial_lfn>22);
    mapped_wm_metric_mask=find(mapped_wm_lfn>19);
    metric_mask=union(pial_metric_mask,mapped_wm_metric_mask);

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
        depth_diff=depth_diff./mapped_wm_depth.*100.0;
    end

     % Add depth and LFN to mask
    if params.perc && params.perc_limit>0
        depth_mask=intersect(mask, find(abs(depth_diff)<=params.perc_limit));    
        all_depth_pial_wm_diff_data(end+1:end+length(depth_mask))=pial_wm_diff_cv(depth_mask);
        all_depth_diff(end+1:end+length(depth_mask))=depth_diff(depth_mask);
        all_depth_contrast(end+1:end+length(depth_mask))=c_idx.*ones(1,length(depth_mask));

        lfn_mask=intersect(mask, find(abs(lfn_diff)<=params.perc_limit));    
        all_lfn_pial_wm_diff_data(end+1:end+length(lfn_mask))=pial_wm_diff_cv(lfn_mask);
        all_lfn_diff(end+1:end+length(lfn_mask))=lfn_diff(lfn_mask);
        all_lfn_contrast(end+1:end+length(lfn_mask))=c_idx.*ones(1,length(lfn_mask));
    else
        all_depth_pial_wm_diff_data(end+1:end+length(mask))=pial_wm_diff_cv(mask);
        all_depth_diff(end+1:end+length(mask))=depth_diff(mask);
        all_depth_contrast(end+1:end+length(mask))=c_idx.*ones(1,length(mask));

        all_lfn_pial_wm_diff_data(end+1:end+length(mask))=pial_wm_diff_cv(mask);
        all_lfn_diff(end+1:end+length(mask))=lfn_diff(mask);
        all_lfn_contrast(end+1:end+length(mask))=c_idx.*ones(1,length(mask));
    end    
end

depth_intercepts=[];
fig=figure();
h=scatterhist(all_depth_diff,all_depth_pial_wm_diff_data,'Group',all_depth_contrast,...
    'Kernel','on','Location','SouthEast','Direction','out','Marker','.',...
    'MarkerSize',5,'Color',colors,'LineStyle',{'-'},'LineWidth',1);
axes(h(1));
hold all;
contrast_names={};
for c_idx=1:length(contrast_order)
    contrast=contrasts(contrast_order(c_idx));
    contrast_names{c_idx}=contrast.comparison_name;
    
    depth_diff=all_depth_diff(find(all_depth_contrast==c_idx));
    pial_wm_diff_data=all_depth_pial_wm_diff_data(find(all_depth_contrast==c_idx));
    
    pPoly = polyfit(depth_diff, pial_wm_diff_data,1); % Linear fit of xdata vs ydata
    depth_intercepts(c_idx)=pPoly(2);
    depth_means(c_idx)=mean(pial_wm_diff_data);
    
    linePointsX = [min(depth_diff) max(depth_diff)]; % find left and right x values
    linePointsY = polyval(pPoly,linePointsX); 
    plot (linePointsX,linePointsY,styles{c_idx},'Color',colors(c_idx,:),'LineWidth',2);
end
legend(contrast_names);
main_xl=xlim();
main_yl=ylim();
axes(h(2));
hold all;
yl=ylim();
for c_idx=1:length(contrast_order)
    idx=find(all_depth_contrast==c_idx);
    contrast_mean=median(all_depth_diff(idx));
    plot([contrast_mean contrast_mean],yl,'--','Color',colors(c_idx,:),'LineWidth',2);
end
xlim(main_xl);
axes(h(3));
hold all;
yl=ylim();
for c_idx=1:length(contrast_order)
    idx=find(all_depth_contrast==c_idx);
    contrast_mean=median(all_depth_pial_wm_diff_data(idx));
    plot([contrast_mean contrast_mean],yl,'--','Color',colors(c_idx,:),'LineWidth',2);
end
xlim(main_yl);
xlabel('Pial depth - WM depth','FontSize',14,'FontName','Arial');
ylabel('Effect size(|Pial Power D| - |WM Power D|)','FontSize',14,'FontName','Arial');
set(gca,'FontName','Arial');

lfn_intercepts=[];
fig=figure();
h=scatterhist(all_lfn_diff,all_lfn_pial_wm_diff_data,'Group',all_lfn_contrast,...
    'Kernel','on','Location','SouthEast','Direction','out','Marker','.',...
    'MarkerSize',5,'Color',colors,'LineStyle',{'-'},'LineWidth',1);
axes(h(1));
hold all;
contrast_names={};
for c_idx=1:length(contrast_order)
    contrast=contrasts(contrast_order(c_idx));
    contrast_names{c_idx}=contrast.comparison_name;
    
    lfn_diff=all_lfn_diff(find(all_lfn_contrast==c_idx));
    pial_wm_diff_data=all_lfn_pial_wm_diff_data(find(all_lfn_contrast==c_idx));
    
    pPoly = polyfit(lfn_diff, pial_wm_diff_data,1); % Linear fit of xdata vs ydata
    lfn_intercepts(c_idx)=pPoly(2);
    lfn_means(c_idx)=mean(pial_wm_diff_data);
    linePointsX = [min(lfn_diff) max(lfn_diff)]; % find left and right x values
    linePointsY = polyval(pPoly,linePointsX); 
    plot (linePointsX,linePointsY,styles{c_idx},'Color',colors(c_idx,:),'LineWidth',2);
end
legend(contrast_names);
main_xl=xlim();
main_yl=ylim();
axes(h(2));
hold all;
yl=ylim();
for c_idx=1:length(contrast_order)
    idx=find(all_lfn_contrast==c_idx);
    contrast_mean=median(all_lfn_diff(idx));
    plot([contrast_mean contrast_mean],yl,'--','Color',colors(c_idx,:),'LineWidth',2);
end
xlim(main_xl);
axes(h(3));
hold all;
yl=ylim();
for c_idx=1:length(contrast_order)
    idx=find(all_lfn_contrast==c_idx);
    contrast_mean=median(all_lfn_pial_wm_diff_data(idx));
    plot([contrast_mean contrast_mean],yl,'--','Color',colors(c_idx,:),'LineWidth',2);
end
xlim(main_yl);
xlabel('Pial LFN - WM LFN','FontSize',14,'FontName','Arial');
ylabel('Effect size(|Pial Power D| - |WM Power D|)','FontSize',14,'FontName','Arial');
set(gca,'FontName','Arial');

rmpath('D:\meg_laminar\layer_comparison');