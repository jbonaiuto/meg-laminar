% Comparison wois - matrix, nx2, n=number of comparisons
% Baseline wois - matrix, nx2, n=number of comparisons
function run_session_foi_layer_comparison(subj_info, session_num,...
    contrast, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'inv_type', 'EBB',...
    'patch_size',0.4, 'surf_dir', '', 'mri_dir', '', 'invert',true,...
    'extract', true, 'compare', true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end
if length(params.surf_dir)==0
    params.surf_dir=fullfile(params.data_dir,'surf');
end

spm('defaults','eeg');
if params.invert
    invert_grey(subj_info, session_num, contrast,...
        'data_dir', params.data_dir, 'inv_type', params.inv_type,...
        'patch_size',params.patch_size, 'surf_dir', params.surf_dir,...
        'mri_dir', params.mri_dir);
end

grey_coreg_dir=fullfile(params.data_dir,'analysis', subj_info.subj_id, num2str(session_num), 'grey_coreg');
foi_dir=fullfile(grey_coreg_dir, params.inv_type,...
    ['p' num2str(params.patch_size)], contrast.zero_event,...
    ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);

if params.extract
    extract_inversion_source(subj_info, session_num, contrast,...
        foi_dir, 'data_dir', params.data_dir, 'inv_type', params.inv_type,...
        'patch_size',params.patch_size,'surf_dir',params.surf_dir);
end

if params.compare
    compare_session_layers(subj_info, session_num, contrast, ...
        foi_dir, 'data_dir', params.data_dir, 'inv_type',params.inv_type,...
        'patch_size', params.patch_size, 'surf_dir',params.surf_dir);
end
