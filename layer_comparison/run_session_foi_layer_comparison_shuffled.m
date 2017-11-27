
function run_session_foi_layer_comparison_shuffled(subj_info, session_num,...
    contrast, shuff_idx, varargin)

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
    invert_grey_shuffled(subj_info, session_num, contrast,...
        shuff_idx, 'data_dir', params.data_dir, 'inv_type', params.inv_type,...
        'patch_size',params.patch_size, 'surf_dir', params.surf_dir,...
        'mri_dir', params.mri_dir);
end

grey_coreg_dir=fullfile(params.data_dir,'analysis', subj_info.subj_id, num2str(session_num), 'grey_coreg');
foi_dir=fullfile(grey_coreg_dir, params.inv_type,...
    ['p' num2str(params.patch_size)], contrast.zero_event,...
    ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))],...
    'shuffled', num2str(shuff_idx));

if params.extract
    extract_inversion_source(subj_info, session_num, contrast,...
        foi_dir, 'data_dir', params.data_dir, 'inv_type', params.inv_type,...
        'patch_size',params.patch_size,'surf_dir',params.surf_dir);
end

if params.compare
    compare_session_layers(subj_info, session_num, contrast, foi_dir,...
        'data_dir', params.data_dir, 'inv_type',params.inv_type,...
        'patch_size', params.patch_size, 'surf_dir',params.surf_dir);
end
rmdir(fullfile(foi_dir, sprintf('t%d_%d', contrast.comparison_woi(1), contrast.comparison_woi(2))),'s');
rmdir(fullfile(foi_dir, sprintf('t%d_%d', contrast.baseline_woi(1), contrast.baseline_woi(2))),'s');
delete(fullfile(foi_dir,sprintf('br%s_%d.dat',subj_info.subj_id,session_num)));
delete(fullfile(foi_dir,sprintf('br%s_%d.mat',subj_info.subj_id,session_num)));
delete(fullfile(foi_dir,sprintf('br%s_%d.surf.gii',subj_info.subj_id,session_num)));
delete(fullfile(foi_dir,sprintf('SPMgainmatrix_br%s_%d_1.mat',subj_info.subj_id,session_num)));
