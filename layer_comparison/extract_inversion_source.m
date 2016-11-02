function extract_inversion_source(subj_info, session_num, foi, woi1, woi2, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'patch_size',0.4);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

grey_coreg_dir=fullfile(params.data_dir,'analysis', subj_info.subj_id, num2str(session_num), 'grey_coreg');
foi_dir=fullfile(grey_coreg_dir, ['p' num2str(params.patch_size)], ['f' num2str(foi(1)) '_' num2str(foi(2))]);

coreg_file_name=fullfile(foi_dir, sprintf('r%s_%d.mat', subj_info.subj_id, session_num));

spm('defaults', 'EEG');
spm_jobman('initcfg'); 

spm_jobman('initcfg'); 
clear jobs

jobfile = 'inversion_results_job.m';
inputs = {};
inputs{1,1}={coreg_file_name};
inputs{2,1}=woi1;
inputs{3,1}=foi;
inputs{4,1}={coreg_file_name};
inputs{5,1}=woi2;
inputs{6,1}=foi;
spm_jobman('initcfg');
spm_jobman('run', jobfile, inputs{:});    

% Move files to subdirectories
woi1_dir=fullfile(foi_dir, ['t' num2str(woi1(1)) '_' num2str(woi1(2))]);
if exist(woi1_dir,'dir')~=7
    mkdir(woi1_dir);
end
woi2_dir=fullfile(foi_dir, ['t' num2str(woi2(1)) '_' num2str(woi2(2))]);
if exist(woi2_dir,'dir')~=7
    mkdir(woi2_dir);
end
movefile(fullfile(foi_dir, sprintf('r%s_%d_1_t%d_%d_f%d_%d_*', subj_info.subj_id, session_num, woi1(1), woi1(2), foi(1), foi(2))), woi1_dir);
movefile(fullfile(foi_dir, sprintf('r%s_%d_1_t%d_%d_f%d_%d_*', subj_info.subj_id, session_num, woi2(1), woi2(2), foi(1), foi(2))), woi2_dir);

% Split pial and grey sources
split_inversion_results(subj_info, grey_coreg_dir, foi, woi1, 'patch_size', params.patch_size);
split_inversion_results(subj_info, grey_coreg_dir, foi, woi2, 'patch_size', params.patch_size);


