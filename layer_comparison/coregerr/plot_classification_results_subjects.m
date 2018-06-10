function plot_classification_results_subjects(subjects, contrast, varargin)

% Parse inputs
defaults = struct('data_dir','d:/meg_laminar/derivatives/spm12',...
    'surf_dir', 'D:/meg_laminar/derivatives/freesurfer','inv_type','EBB',...
    'patch_size',0.4,'recompute_roi',false,'iterations',10,'shift_magnitude',10);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end
addpath('D:\meg_laminar\layer_comparison');
spm('defaults','eeg');

subject_whole_brain_tvals=zeros(length(subjects), params.iterations);
subject_func_roi_tvals=zeros(length(subjects), params.iterations);
subject_anatfunc_roi_tvals=zeros(length(subjects), params.iterations);

thresh_type='lower';
switch contrast.comparison_name
    case 'dots_beta_erd'
        thresh_type='upper';
    case 'dots_alpha'
        thresh_type='upper';
end

subj_dofs=zeros(length(subjects), params.iterations);

for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    surf_dir=fullfile(params.surf_dir, subj_info.subj_id);
    orig_white_mesh=fullfile(surf_dir,'white.hires.deformed.surf.gii');
    white_mesh=fullfile(surf_dir,'ds_white.hires.deformed.surf.gii');
    white_inflated=fullfile(surf_dir,'ds_white.hires.deformed_inflated.surf.gii');
    orig_pial_mesh=fullfile(surf_dir,'pial.hires.deformed.surf.gii');
    pial_mesh=fullfile(surf_dir,'ds_pial.hires.deformed.surf.gii');
    pial_inflated=fullfile(surf_dir,'ds_pial.hires.deformed_inflated.surf.gii');
    pial_white_map=map_pial_to_white(white_mesh, pial_mesh, 'mapType', 'link',...
        'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh);
    white_pial_map=map_white_to_pial(white_mesh, pial_mesh, 'mapType', 'link',...
        'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh);
    pial_hemisphere_map=get_hemisphere_map(pial_mesh, orig_pial_mesh);
    white_hemisphere_map=get_hemisphere_map(white_mesh, orig_white_mesh);
    
    foi_dir=fullfile(params.data_dir, subj_info.subj_id,...
            sprintf('ses-%02d',subj_info.sessions(1)), 'grey_coreg', params.inv_type,....
            ['p' num2str(params.patch_size)], contrast.zero_event,...
            ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
    lfn_filename=fullfile(foi_dir, sprintf('br%s_%d.mat',subj_info.subj_id, subj_info.sessions(1)));    
    
    orig_foi_dir=fullfile(params.data_dir, subj_info.subj_id,...
                'grey_coreg', params.inv_type,....
                ['p' num2str(params.patch_size)], contrast.zero_event,...
                ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
    
    for idx=1:params.iterations
        foi_dir=fullfile(params.data_dir, subj_info.subj_id,...
                'grey_coreg', params.inv_type,....
                ['p' num2str(params.patch_size)], contrast.zero_event,...
                ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))],...
                'coregerr', num2str(params.shift_magnitude), num2str(idx));
        % Get mean pial-wm diff
        pial_wm_diff=gifti(fullfile(foi_dir,['pial-white.' contrast.comparison_name '.diff.gii']));

        subj_dofs(subj_idx,idx)=size(pial_wm_diff.cdata(:,:),2)-1;
        
        % whole brain                
        [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, contrast.comparison_name, ...
            thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
            pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',0,...
            'type','mean', 'region', '', 'hemisphere', '',...
            'pial_hemisphere_map', pial_hemisphere_map,...
            'white_hemisphere_map', white_hemisphere_map, 'recompute', params.recompute_roi);            
        pial_wm_roi_diff=mean(pial_wm_diff.cdata(mask,:));
        % Perform ROI t-stat
        [tstat,p]=ttest_corrected(pial_wm_roi_diff','correction',25*var(pial_wm_roi_diff));
        disp(sprintf('%s, %d, whole=%.3f',subj_info.subj_id,idx,tstat));
        subject_whole_brain_tvals(subj_idx,idx)=tstat;

        % func ROI
        thresh=80;
        [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, contrast.comparison_name, ...
            thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
            pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',thresh,...
            'type','mean', 'region', '', 'hemisphere', '',...
            'pial_hemisphere_map', pial_hemisphere_map,...
            'white_hemisphere_map', white_hemisphere_map, 'recompute', params.recompute_roi);            
        pial_wm_roi_diff=mean(pial_wm_diff.cdata(mask,:));
        % Perform ROI t-stat
        [tstat,p]=ttest_corrected(pial_wm_roi_diff','correction',25*var(pial_wm_roi_diff));
        disp(sprintf('%s, %d, func=%.3f',subj_info.subj_id,idx,tstat));
        subject_func_roi_tvals(subj_idx,idx)=tstat;

        % anat-func ROI
        [pial_mask,wm_mask,mask]=compute_roi(subj_info, orig_foi_dir, contrast.comparison_name, ...
            thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
            pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',thresh,...
            'type','mean', 'region', contrast.region, 'hemisphere', contrast.hemisphere,...
            'pial_hemisphere_map', pial_hemisphere_map,...
            'white_hemisphere_map', white_hemisphere_map, 'recompute', params.recompute_roi);            
        pial_wm_roi_diff=mean(pial_wm_diff.cdata(mask,:));
        % Perform ROI t-stat
        [tstat,p]=ttest_corrected(pial_wm_roi_diff','correction',25*var(pial_wm_roi_diff));
        disp(sprintf('%s, %d, anat-func=%.3f',subj_info.subj_id,idx,tstat));
        subject_anatfunc_roi_tvals(subj_idx,idx)=tstat;    
    end
end

fig=figure('Position',[1 1 1250 600],'PaperUnits','points',...
    'PaperPosition',[1 1 650 300],'PaperPositionMode','manual');
hold on;

subj_ids={};
alpha=1.0-(0.05/2);
t_thresh=tinv(alpha, mean(subj_dofs(:)));
plot([1-.5 length(subjects)+.75],[t_thresh t_thresh],'k--','LineWidth',2);
plot([1-.5 length(subjects)+.75],[-t_thresh -t_thresh],'k--','LineWidth',2);

for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    subj_ids{subj_idx}=num2str(subj_idx);
    
    for idx=1:params.iterations   
        % Whole brain bar
        tval=subject_whole_brain_tvals(subj_idx,idx);
        plot(subj_idx-.2,tval,'o');

        % Func ROI bar
        tval=subject_func_roi_tvals(subj_idx,idx);
        plot(subj_idx,tval,'*');

        % Anat-Func ROI bar
        tval=subject_anatfunc_roi_tvals(subj_idx,idx);
        plot(subj_idx+.2,tval,'+');
    end
end
set(gca,'XTick',1:length(subjects));
set(gca,'XTickLabel',subj_ids);
xlabel('Participant','Fontsize',24,'Fontname','Arial');
ylabel('Pial-White ROI t-statistic','Fontsize',24,'Fontname','Arial');
xlim([.5 8.5]);
yl=ylim;
ylim([-max(abs(yl)) max(abs(yl))]);
xt=get(gca,'XTick');
set(gca,'FontSize',20);
set(gca,'Fontname','Arial');

rmpath('D:\meg_laminar\layer_comparison');