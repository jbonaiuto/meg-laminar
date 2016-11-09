function all_subjects_sensor_tf(subjects, zero_event, varargin)

defaults = struct('data_dir', '/data/pred_coding');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults', 'EEG');
spm_jobman('initcfg');

cd D:\pred_coding\src\matlab\analysis\sensor_tf
spm_jobman('initcfg');
clear jobs;    
jobfile = 'second_level_job.m';
% Scalp x freq
inputs={};
inputs{1,1}={fullfile(params.data_dir, 'analysis', ['scalp_freq_rtf_rc' zero_event '_Tafdf_positive'])};  
scans={};
for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    scans{end+1,1}=fullfile(params.data_dir, 'analysis', subj_info.subj_id, ['scalp_freq_rtf_rc' zero_event '_Tafdf'], 'con_0005.nii,1');
end
inputs{2,1}=scans;
inputs{3,1}=3;
spm_jobman('run', jobfile, inputs{:});

cd D:\pred_coding\src\matlab\analysis\sensor_tf
spm_jobman('initcfg');
clear jobs;    
jobfile = 'second_level_job.m';
% Scalp x freq
inputs={};
inputs{1,1}={fullfile(params.data_dir, 'analysis', ['scalp_freq_rtf_rc' zero_event '_Tafdf_negative'])};  
scans={};
for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    scans{end+1,1}=fullfile(params.data_dir, 'analysis', subj_info.subj_id, ['scalp_freq_rtf_rc' zero_event '_Tafdf'], 'con_0011.nii,1');
end
inputs{2,1}=scans;
inputs{3,1}=3;
spm_jobman('run', jobfile, inputs{:});

cd D:\pred_coding\src\matlab\analysis\sensor_tf
spm_jobman('initcfg');
clear jobs;    
jobfile = 'second_level_job.m';
% Scalp x time
inputs={};
inputs{1,1}={fullfile(params.data_dir, 'analysis', ['scalp_time_rtf_rc' zero_event '_Tafdf_positive'])};  
scans={};
for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    scans{end+1,1}=fullfile(params.data_dir, 'analysis', subj_info.subj_id, ['scalp_time_rtf_rc' zero_event '_Tafdf'], 'con_0005.nii,1');
end
inputs{2,1}=scans;
inputs{3,1}=2;
spm_jobman('run', jobfile, inputs{:});

cd D:\pred_coding\src\matlab\analysis\sensor_tf
spm_jobman('initcfg');
clear jobs;    
jobfile = 'second_level_job.m';
% Scalp x time
inputs={};
inputs{1,1}={fullfile(params.data_dir, 'analysis', ['scalp_time_rtf_rc' zero_event '_Tafdf_negative'])};  
scans={};
for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    scans{end+1,1}=fullfile(params.data_dir, 'analysis', subj_info.subj_id, ['scalp_time_rtf_rc' zero_event '_Tafdf'], 'con_0011.nii,1');
end
inputs{2,1}=scans;
inputs{3,1}=2;
spm_jobman('run', jobfile, inputs{:});

cd D:\pred_coding\src\matlab\analysis\sensor_tf
spm_jobman('initcfg');
clear jobs;    
jobfile = 'second_level_job.m';
% Scalp x time
inputs={};
inputs{1,1}={fullfile(params.data_dir, 'analysis', ['time_freq_rtf_rc' zero_event '_Tafdf_positive'])};  
scans={};
for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    scans{end+1,1}=fullfile(params.data_dir, 'analysis', subj_info.subj_id, ['time_freq_rtf_rc' zero_event '_Tafdf'], 'con_0005.nii,1');
end
inputs{2,1}=scans;
inputs{3,1}=4;
spm_jobman('run', jobfile, inputs{:});

cd D:\pred_coding\src\matlab\analysis\sensor_tf
spm_jobman('initcfg');
clear jobs;    
jobfile = 'second_level_job.m';
% Scalp x time
inputs={};
inputs{1,1}={fullfile(params.data_dir, 'analysis', ['time_freq_rtf_rc' zero_event '_Tafdf_negative'])};  
scans={};
for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    scans{end+1,1}=fullfile(params.data_dir, 'analysis', subj_info.subj_id, ['time_freq_rtf_rc' zero_event '_Tafdf'], 'con_0011.nii,1');
end
inputs{2,1}=scans;
inputs{3,1}=4;
spm_jobman('run', jobfile, inputs{:});
cd D:\pred_coding\src\matlab\analysis\sensor_tf