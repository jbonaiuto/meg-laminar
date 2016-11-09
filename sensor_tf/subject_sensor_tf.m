function subject_sensor_tf(subj_info, zero_event, varargin)

defaults = struct('data_dir', '/data/pred_coding');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults', 'EEG');
spm_jobman('initcfg');

conditions={'congruent-low','congruent-med','congruent-high',...
    'incongruent-low','incongruent-med','incongruent-high'};

cd D:\pred_coding\src\matlab\analysis\sensor_tf
spm_jobman('initcfg');
clear jobs;    
jobfile = 'factorial_job.m';
% Scalp x freq
inputs={};
inputs{1,1}={fullfile(params.data_dir, 'analysis', ...
        subj_info.subj_id, ['scalp_freq_rtf_rc' zero_event '_Tafdf'])};    
for cond_idx=1:length(conditions)
    level_scans={};
    for session_num=1:length(subj_info.sessions)
        img_path=fullfile(params.data_dir, 'analysis', ...
            subj_info.subj_id, num2str(session_num), ...
            ['scalp_freq_rtf_rc' zero_event '_Tafdf' num2str(session_num)],...
            ['scondition_' conditions{cond_idx} '.nii']);
        V=spm_vol(img_path);
        for i=1:length(V)
            level_scans{end+1,1}=[img_path ',' num2str(i)];
        end
    end
    inputs{1+cond_idx,1}=level_scans;
end
inputs{8,1}=3;
inputs{9,1}=3;
inputs{10,1}=3;
spm_jobman('run', jobfile, inputs{:});

% Scalp x time
cd D:\pred_coding\src\matlab\analysis\sensor_tf
spm_jobman('initcfg');
clear jobs;    
inputs={};
inputs{1,1}={fullfile(params.data_dir, 'analysis', ...
        subj_info.subj_id, ['scalp_time_rtf_rc' zero_event '_Tafdf'])};    
for cond_idx=1:length(conditions)
    level_scans={};
    for session_num=1:length(subj_info.sessions)
        img_path=fullfile(params.data_dir, 'analysis', ...
            subj_info.subj_id, num2str(session_num), ...
            ['scalp_time_rtf_rc' zero_event '_Tafdf' num2str(session_num)],...
            ['scondition_' conditions{cond_idx} '.nii']);
        V=spm_vol(img_path);
        for i=1:length(V)
            level_scans{end+1,1}=[img_path ',' num2str(i)];
        end
    end
    inputs{1+cond_idx,1}=level_scans;
end
inputs{8,1}=2;
inputs{9,1}=2;
inputs{10,1}=2;
spm_jobman('run', jobfile, inputs{:});

% Time x frequency
cd D:\pred_coding\src\matlab\analysis\sensor_tf
spm_jobman('initcfg');
clear jobs;    
inputs={};
inputs{1,1}={fullfile(params.data_dir, 'analysis', ...
        subj_info.subj_id, ['time_freq_rtf_rc' zero_event '_Tafdf'])};    
for cond_idx=1:length(conditions)
    level_scans={};
    for session_num=1:length(subj_info.sessions)
        img_path=fullfile(params.data_dir, 'analysis', ...
            subj_info.subj_id, num2str(session_num), ...
            ['time_freq_rtf_rc' zero_event '_Tafdf' num2str(session_num)],...
            ['scondition_' conditions{cond_idx} '.nii']);
        V=spm_vol(img_path);
        for i=1:length(V)
            level_scans{end+1,1}=[img_path ',' num2str(i)];
        end        
    end
    inputs{1+cond_idx,1}=level_scans;
end
inputs{8,1}=4;
inputs{9,1}=4;
inputs{10,1}=4;
spm_jobman('run', jobfile, inputs{:});

cd D:\pred_coding\src\matlab\analysis\sensor_tf
