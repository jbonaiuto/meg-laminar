function extract_inversion_source(subj_info, session_num, foi, woi, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'inv_type', 'EBB', 'patch_size',0.4);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

grey_coreg_dir=fullfile(params.data_dir,'analysis', subj_info.subj_id, num2str(session_num), 'grey_coreg');
foi_dir=fullfile(grey_coreg_dir, params.inv_type, ['p' num2str(params.patch_size)], ['f' num2str(foi(1)) '_' num2str(foi(2))]);

coreg_file_name=fullfile(foi_dir, sprintf('r%s_%d.mat', subj_info.subj_id, session_num));

spm('defaults', 'EEG');

% Extract source at woi1 and woi2
spm_jobman('initcfg'); 
clear jobs
matlabbatch={};
batch_idx=1;

matlabbatch{batch_idx}.spm.meeg.source.results.D = {coreg_file_name};
matlabbatch{batch_idx}.spm.meeg.source.results.val = 1;
matlabbatch{batch_idx}.spm.meeg.source.results.woi = woi;
matlabbatch{batch_idx}.spm.meeg.source.results.foi = foi;
matlabbatch{batch_idx}.spm.meeg.source.results.ctype = 'trials';
matlabbatch{batch_idx}.spm.meeg.source.results.space = 0;
matlabbatch{batch_idx}.spm.meeg.source.results.format = 'mesh';
matlabbatch{batch_idx}.spm.meeg.source.results.smoothing = 8;
batch_idx=batch_idx+1;

spm_jobman('run',matlabbatch);

% Move files to subdirectories
woi_dir=fullfile(foi_dir, ['t' num2str(woi(1)) '_' num2str(woi(2))]);
if exist(woi_dir,'dir')~=7
    mkdir(woi_dir);
end
movefile(fullfile(foi_dir, ['r' subj_info.subj_id '_' num2str(session_num) '_1_t' num2str(woi(1)) '_' num2str(woi(2)) '_f' num2str(foi(1)) '_' num2str(foi(2)) '_*']), woi_dir);

% Split pial and grey sources
split_inversion_results(subj_info, grey_coreg_dir, foi, woi, 'patch_size', params.patch_size, 'inv_type', params.inv_type);




