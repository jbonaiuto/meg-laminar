function [pial_roi_trials,wm_roi_trials,pial_roi_bc_trials,...
    wm_roi_bc_trials,times,freqs]=compute_condition_power(subj_info,...
    session_num, contrast, varargin)

% Parse inputs
defaults = struct('data_dir','d:/pred_coding/derivatives/spm12','surf_dir', 'D:/pred_coding/derivatives/freesurfer',...
    'inv_type','EBB','patch_size',0.4,'thresh_percentile',80,'roi_type','mean','recompute',false,...
    'recompute_roi',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

freqs=[contrast.foi(1)-5:contrast.foi(2)+5];

grey_coreg_dir=fullfile(params.data_dir,subj_info.subj_id,...
    sprintf('ses-%02d', session_num), 'grey_coreg');
foi_dir=fullfile(grey_coreg_dir, params.inv_type, ...
    ['p' num2str(params.patch_size)], contrast.zero_event,...
    ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);

pial_data_filename=fullfile(foi_dir, sprintf('pial.%s.roi_rtf.mat', contrast.comparison_name));
wm_data_filename=fullfile(foi_dir, sprintf('white.%s.roi_rtf.mat', contrast.comparison_name));

if (exist(pial_data_filename,'file')==2 && exist(wm_data_filename,'file')==2) && ~params.recompute
    load(pial_data_filename);
    pial_roi_bc_trials=results.pial_roi_bc_trials;
    times=results.times;
    freqs=results.freqs;
    load(wm_data_filename);
    wm_roi_bc_trials=results.wm_roi_bc_trials;
    
    pial_data_filename=fullfile(foi_dir, sprintf('pial.%s.roi_tf.mat', contrast.comparison_name));
    load(pial_data_filename);
    pial_roi_trials=results.pial_roi_trials;
    wm_data_filename=fullfile(foi_dir, sprintf('white.%s.roi_tf.mat', contrast.comparison_name));
    load(wm_data_filename);
    wm_roi_trials=results.wm_roi_trials;
else
    addpath('../layer_comparison');

    thresh_type='lower';
    switch contrast.direction
        case 'negative'
            thresh_type='upper';
        case 'positive'
            thresh_type='lower';
    end

    spm('defaults','eeg');

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
            sprintf('ses-%02d', session_num), 'grey_coreg', params.inv_type,....
            ['p' num2str(params.patch_size)], contrast.zero_event,...
            ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
    lfn_filename=fullfile(foi_dir, sprintf('br%s_%d.mat',subj_info.subj_id, session_num));

    grey_coreg_dir=fullfile(params.data_dir, subj_info.subj_id,...
        sprintf('ses-%02d', session_num), 'grey_coreg');
    foi_dir=fullfile(grey_coreg_dir, params.inv_type,...
        ['p' num2str(params.patch_size)], contrast.zero_event,...
        ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);

    % Get ROI
    [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, contrast.comparison_name, ...
                thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
                pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',params.thresh_percentile,...
                'type',params.roi_type, 'region', contrast.region, 'hemisphere', contrast.hemisphere,...
                'pial_hemisphere_map', pial_hemisphere_map,...
                'white_hemisphere_map', white_hemisphere_map, 'recompute', params.recompute_roi);
    mapped_wm_mask=pial_white_map(mask);
    
    % Load coreg file
    coreg_file_name=fullfile(foi_dir, sprintf('br%s_%d.mat', subj_info.subj_id, session_num));
    D=spm_eeg_load(coreg_file_name);
    goodchans=D.indchantype('MEGGRAD','good');
    M=D.inv{1}.inverse.M;
    U=D.inv{1}.inverse.U{1};
    nvertices=size(M,1)/2;
    It   = D.inv{1}.inverse.It;
    times=D.inv{1}.inverse.pst;
    Dgood=squeeze(D(goodchans,It,:));
    ntrials=size(Dgood,3);
    pial_MU=M(mask+nvertices,:)*U;
    wm_MU=M(mapped_wm_mask,:)*U;

    % Wavelet for each vertex in ROI for each trial
    S.frequencies=freqs;
    pial_roi_trials=zeros(length(S.frequencies),length(times),ntrials);
    wm_roi_trials=zeros(length(S.frequencies),length(times),ntrials);

    
    baseline_idx=intersect(find(times>=contrast.baseline_woi(1)),find(times<=contrast.baseline_woi(2)));

    for t=1:ntrials
        t
        pial_MUtrial=pial_MU*squeeze(Dgood(:,:,t));
        wm_MUtrial=wm_MU*squeeze(Dgood(:,:,t));

        pial_trial=spm_eeg_specest_morlet(S,pial_MUtrial,times.*1e-3);
        pial_trial_tf=pial_trial.fourier.*conj(pial_trial.fourier);
        pial_roi_trials(:,:,t)=squeeze(mean(pial_trial_tf));


        wm_trial=spm_eeg_specest_morlet(S,wm_MUtrial,times.*1e-3);
        wm_trial_tf=wm_trial.fourier.*conj(wm_trial.fourier);
        wm_roi_trials(:,:,t)=squeeze(mean(wm_trial_tf));
    end

    pial_data_filename=fullfile(foi_dir, sprintf('pial.%s.roi_tf.mat', contrast.comparison_name));
    wm_data_filename=fullfile(foi_dir, sprintf('white.%s.roi_tf.mat', contrast.comparison_name));

    results=[];
    results.pial_roi_trials=pial_roi_trials;
    results.times=times;
    results.freqs=freqs;
    save(pial_data_filename, 'results');
    
    results=[];
    results.wm_roi_trials=wm_roi_trials;
    results.times=times;
    results.freqs=freqs;
    save(wm_data_filename, 'results');
    
    pial_roi_bc_trials=zeros(length(S.frequencies),length(times),ntrials);
    wm_roi_bc_trials=zeros(length(S.frequencies),length(times),ntrials);
    
    for t=1:ntrials
        pial_trial=squeeze(pial_roi_trials(:,:,t));
        pial_trial_baseline = spm_robust_average(pial_trial(:,baseline_idx),2);
        pial_trial_baseline = repmat(pial_trial_baseline,[1 length(times)]);
        pial_roi_bc_trials(:,:,t) = (pial_trial - pial_trial_baseline)./pial_trial_baseline.*100;
        
        wm_trial=squeeze(wm_roi_trials(:,:,t));
        wm_trial_baseline = spm_robust_average(wm_trial(:,baseline_idx),2);
        wm_trial_baseline=repmat(wm_trial_baseline,[1 length(times)]);
        wm_roi_bc_trials(:,:,t) = (wm_trial - wm_trial_baseline)./wm_trial_baseline.*100;
    end
    
    pial_data_filename=fullfile(foi_dir, sprintf('pial.%s.roi_rtf.mat', contrast.comparison_name));
    wm_data_filename=fullfile(foi_dir, sprintf('white.%s.roi_rtf.mat', contrast.comparison_name));

    results=[];
    results.pial_roi_bc_trials=pial_roi_bc_trials;
    results.times=times;
    results.freqs=freqs;
    save(pial_data_filename, 'results');
    
    results=[];
    results.wm_roi_bc_trials=wm_roi_bc_trials;
    results.times=times;
    results.freqs=freqs;
    save(wm_data_filename, 'results');

end
