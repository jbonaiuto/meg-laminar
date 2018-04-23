function run_all_session(subj_info, session_num, contrasts, varargin)

defaults = struct('data_dir', 'd:/pred_coding', 'inv_type', 'EBB',...
    'patch_size',0.4, 'surf_dir', '', 'mri_dir', '',... 
    'invert',true, 'extract', true, 'compare', true, 'shuffle_iterations',10);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end
if length(params.mri_dir)==0
    params.mri_dir=fullfile(params.data_dir,'mri');
end
if length(params.surf_dir)==0
    params.surf_dir=fullfile(params.data_dir,'surf');
end

for shuff_idx=1:params.shuffle_iterations
    % Redo shuffling
    clear spm_eeg_invert_classic;
    for i=1:length(contrasts)
        contrast=contrasts(i);
        run_session_foi_layer_comparison(subj_info, session_num,...
            contrast, shuff_idx, 'data_dir', params.data_dir, 'inv_type', params.inv_type,...
            'patch_size',params.patch_size, 'surf_dir', params.surf_dir,...
            'mri_dir', params.mri_dir, 'invert',params.invert,'extract', params.extract,...
            'compare', params.compare);
    end        
end