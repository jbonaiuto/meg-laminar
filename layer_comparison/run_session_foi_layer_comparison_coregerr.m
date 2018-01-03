function run_session_foi_layer_comparison_coregerr(subj_info, session_num,...
    contrast, idx, varargin)

% Parse inputs
defaults = struct('data_dir', 'd:/pred_coding/derivatives/spm12', 'inv_type', 'EBB',...
    'patch_size',0.4, 'surf_dir', 'd:/pred_coding/derivatives/freesurfer', 'mri_dir', 'd:/pred_coding', 'invert',true,...
    'extract', true, 'compare', true, 'shift_magnitude', 10);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

if params.invert
    invert_grey_coregerr(subj_info, session_num, contrast,...
        idx, 'data_dir', params.data_dir, 'inv_type', params.inv_type,...
        'patch_size',params.patch_size, 'surf_dir', params.surf_dir,...
        'mri_dir', params.mri_dir, 'shift_magnitude', params.shift_magnitude);
end

grey_coreg_dir=fullfile(params.data_dir, subj_info.subj_id, sprintf('ses-0%d',session_num), 'grey_coreg');
foi_dir=fullfile(grey_coreg_dir, params.inv_type,...
    ['p' num2str(params.patch_size)], contrast.zero_event,...
    ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))],...);
    'coregerr', num2str(params.shift_magnitude), num2str(idx));

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
rmdir(fullfile(foi_dir, sprintf('t%d_%d', contrast.baseline_woi(1), contrast.baseline_woi(2))),'s');
rmdir(fullfile(foi_dir, sprintf('t%d_%d', contrast.comparison_woi(1), contrast.comparison_woi(2))),'s');
delete(fullfile(foi_dir,sprintf('br%s_%d.dat',subj_info.subj_id,session_num)));
delete(fullfile(foi_dir,sprintf('br%s_%d.mat',subj_info.subj_id,session_num)));
delete(fullfile(foi_dir,sprintf('br%s_%d.surf.gii',subj_info.subj_id,session_num)));
delete(fullfile(foi_dir,sprintf('SPMgainmatrix_br%s_%d_1.mat',subj_info.subj_id,session_num)));
