function run_all_session(subj_info, contrasts, varargin)

defaults = struct('data_dir', 'd:/meg_laminar/derivatives/spm12', 'inv_type', 'EBB',...
    'patch_size',0.4, 'surf_dir', 'd:/meg_laminar/derivatives/freesurfer',...
    'mri_dir', 'd:/meg_laminar', 'invert',true, 'extract', true, 'compare', true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

addpath('D:\meg_laminar\layer_comparison');

for i=1:length(contrasts)
    contrast=contrasts(i);
    for s_idx=1:length(subj_info.sessions)
        run_session_foi_layer_comparison(subj_info, s_idx,...
            contrast, 'data_dir', params.data_dir, 'inv_type', params.inv_type,...
            'patch_size',params.patch_size, 'surf_dir', params.surf_dir,...
            'mri_dir', params.mri_dir, 'invert',params.invert,'extract', params.extract,...
            'compare', params.compare);
    end
end        

rmpath('D:\meg_laminar\layer_comparison');