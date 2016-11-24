function invert_grey(subj_info, session_num, zero_event, foi, woi, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'patch_size',0.4, 'surf_dir', '', 'mri_dir', '');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end
if length(params.surf_dir)==0
    params.surf_dir=fullfile(params.data_dir,'surf');
end
if length(params.mri_dir)==0
    params.mri_dir=fullfile(params.data_dir,'mri');
end

data_dir=fullfile(params.data_dir,'analysis', subj_info.subj_id, num2str(session_num));
data_file_name=fullfile(data_dir, sprintf('rc%s_Tafdf%d.mat', zero_event, session_num));

% Create directory for inversion results
foi_dir=fullfile(data_dir, 'grey_coreg', ['p' num2str(params.patch_size)], ['f' num2str(foi(1)) '_' num2str(foi(2))]);
if exist(foi_dir,'dir')~=7
    mkdir(foi_dir);
end
coreg_file_name=fullfile(foi_dir, sprintf('%s_%d.mat', subj_info.subj_id, session_num));
removed_file_name=fullfile(foi_dir, sprintf('r%s_%d.mat', subj_info.subj_id, session_num));

spm('defaults', 'EEG');
spm_jobman('initcfg'); 

% Copy file to foi_dir
clear jobs
matlabbatch={};
batch_idx=1;
matlabbatch{batch_idx}.spm.meeg.other.copy.D = {data_file_name};
matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = coreg_file_name;
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.meeg.preproc.remove.D = {coreg_file_name};
matlabbatch{batch_idx}.spm.meeg.preproc.remove.prefix = 'r';
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.meeg.preproc.prepare.D = {removed_file_name};
matlabbatch{batch_idx}.spm.meeg.preproc.prepare.task{1}.settype.channels{1}.type = 'EEG';
matlabbatch{batch_idx}.spm.meeg.preproc.prepare.task{1}.settype.newtype = 'Other';
batch_idx=batch_idx+1;
    
spm_jobman('run',matlabbatch);
 
% Relabel trials to all be same condition
load(removed_file_name);
D.condlist={zero_event};
for trial_idx=1:length(D.trials)
    D.trials(trial_idx).label=zero_event;
end
save(removed_file_name,'D');

spm_jobman('initcfg'); 
clear jobs

jobfile = 'invert_grey_job.m';
inputs = {};
inputs{1,1}={removed_file_name};
inputs{2,1}={fullfile(params.mri_dir,[subj_info.subj_id subj_info.birth_date], [subj_info.headcast_t1 ',1'])};
inputs{3,1}={fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed-ds_pial.hires.deformed.surf.gii')};
inputs{4,1}=subj_info.nas;
inputs{5,1}=subj_info.lpa;
inputs{6,1}=subj_info.rpa;
inputs{7,1}=woi;
inputs{8,1}=foi;
inputs{9,1}=params.patch_size;
spm_jobman('initcfg');
spm_jobman('run', jobfile, inputs{:});    
