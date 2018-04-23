function run_all_subject_shuffled(subj_info, contrasts, varargin)

defaults = struct('data_dir', 'd:/pred_coding/derivatives/spm12', 'inv_type', 'EBB',...
    'patch_size',0.4, 'surf_dir', 'd:/pred_coding/derivatives/freesurfer', 'filter_sessions',true, 'shuffle_iterations',10);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

for shuff_idx=1:params.shuffle_iterations
    % Redo shuffling
    clear spm_eeg_invert_classic;
    for i=1:length(contrasts)
        contrast=contrasts(i);
        compare_subject_layers_shuffled(subj_info, contrast, shuff_idx,...
            'data_dir', params.data_dir, 'save_results', true, 'inv_type', params.inv_type, ...
            'patch_size', params.patch_size, 'surf_dir', params.surf_dir,...
            'filter_sessions', params.filter_sessions);
    end        
end
