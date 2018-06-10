function export_classification_results_against_shuffled_subjects(subjects, contrast, varargin)

% Parse inputs
defaults = struct('data_dir','d:/meg_laminar/derivatives/spm12',...
    'surf_dir', 'D:/meg_laminar/derivatives/freesurfer','inv_type','EBB',...
    'patch_size',0.4,'recompute_roi',false, 'iterations', 10);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

fid=fopen(fullfile(params.data_dir, 'against_shuffled',contrast.comparison_name,'subjects_data.csv'),'w');
fprintf(fid, 'Subject,Trial,ROIType,Pial-White\n');

thresh_type='lower';
switch contrast.comparison_name
    case 'dots_beta_erd'
        thresh_type='upper';
    case 'dots_alpha'
        thresh_type='upper';
end
subj_dofs=zeros(1,length(subjects));

for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    surf_dir=fullfile(params.surf_dir,subj_info.subj_id);
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
            sprintf('ses-%02d', subj_info.sessions(1)), 'grey_coreg', params.inv_type,....
            ['p' num2str(params.patch_size)], contrast.zero_event,...
            ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
    lfn_filename=fullfile(foi_dir, sprintf('br%s_%d.mat',subj_info.subj_id, subj_info.sessions(1)));    
    
    pial_wm_diff_shuffled=[];
    for shuf_idx=1:params.iterations
        foi_dir=fullfile(params.data_dir, subj_info.subj_id,...
                'grey_coreg', params.inv_type,....
                ['p' num2str(params.patch_size)], contrast.zero_event,...
                ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))],...
                'shuffled', num2str(shuf_idx));
        % Get mean pial-wm in ROI
        pial_wm_diff=gifti(fullfile(foi_dir,['pial-white.' contrast.comparison_name '.diff.gii']));                
        pial_wm_diff_shuffled(shuf_idx,:,:)=pial_wm_diff.cdata(:,:);
    end
    pial_wm_diff_shuffled=squeeze(mean(pial_wm_diff_shuffled));
    
    foi_dir=fullfile(params.data_dir, subj_info.subj_id,...
            'grey_coreg', params.inv_type,....
            ['p' num2str(params.patch_size)], contrast.zero_event,...
            ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
    % Get mean pial-wm in ROI
    pial_wm_diff=gifti(fullfile(foi_dir,['pial-white.' contrast.comparison_name '.diff.gii']));
            
    subj_dofs(subj_idx)=size(pial_wm_diff.cdata(:,:),2)-1;
    
    % whole brain
    [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, contrast.comparison_name, ...
        thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
        pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',0,...
        'type','mean', 'region', '', 'hemisphere', '',...
        'pial_hemisphere_map', pial_hemisphere_map,...
        'white_hemisphere_map', white_hemisphere_map, 'recompute', params.recompute_roi);            
    pial_wm_roi_diff=mean(pial_wm_diff.cdata(mask,:))-mean(pial_wm_diff_shuffled(mask,:));
    for t_idx=1:length(pial_wm_roi_diff)
        fprintf(fid, sprintf('%d,%d,Global,%.8f\n', subj_idx, t_idx, pial_wm_roi_diff(t_idx)));
    end
    
    % func ROI
    [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, contrast.comparison_name, ...
        thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
        pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',80,...
        'type','mean', 'region', '', 'hemisphere', '',...
        'pial_hemisphere_map', pial_hemisphere_map,...
        'white_hemisphere_map', white_hemisphere_map, 'recompute', params.recompute_roi);            
    pial_wm_roi_diff=mean(pial_wm_diff.cdata(mask,:))-mean(pial_wm_diff_shuffled(mask,:));
    for t_idx=1:length(pial_wm_roi_diff)
        fprintf(fid, sprintf('%d,%d,Func,%.8f\n', subj_idx, t_idx, pial_wm_roi_diff(t_idx)));
    end
    
    % anat-func ROI
    [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, contrast.comparison_name, ...
        thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
        pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',80,...
        'type','mean', 'region', contrast.region, 'hemisphere', contrast.hemisphere,...
        'pial_hemisphere_map', pial_hemisphere_map,...
        'white_hemisphere_map', white_hemisphere_map, 'recompute', params.recompute_roi);            
    pial_wm_roi_diff=mean(pial_wm_diff.cdata(mask,:))-mean(pial_wm_diff_shuffled(mask,:));
    for t_idx=1:length(pial_wm_roi_diff)
        fprintf(fid, sprintf('%d,%d,Anat,%.8f\n', subj_idx, t_idx, pial_wm_roi_diff(t_idx)));
    end
end

fclose(fid);