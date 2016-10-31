% resp woi=[-1500 1500], baseline=[-1500 -1000]
% instr woi=[-3000 1000], baseline=[-3000 -2500]
function session_sensor_tf(subj_info, session_num, zero_event, woi, baseline, varargin)

defaults = struct('data_dir', '/data/pred_coding', 'run_tf', true, ...
    'run_smooth', true, 'run_stats', true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

conditions={'congruent-low','congruent-med','congruent-high',...
    'incongruent-low','incongruent-med','incongruent-high'};

% Run TF decomposition
if params.run_tf
    preprocessed_filename=fullfile(params.data_dir, 'analysis', ...
        subj_info.subj_id, num2str(session_num), ...
        ['rc' zero_event '_Tafdf' num2str(session_num) '.mat']);
    jobfile = 'tf_job.m';
    inputs = {};
    inputs{1,1}={preprocessed_filename};
    inputs{2,1}=baseline;
    inputs{3,1}=woi;
    inputs{4,1}=woi;
    inputs{5,1}=woi;
    spm('defaults', 'EEG');
    spm_jobman('run', jobfile, inputs{:});
end

% Smooth resulting images
if params.run_smooth
    smooth_files={};
    for cond_idx=1:length(conditions)
        img_path=fullfile(params.data_dir, 'analysis', ...
            subj_info.subj_id, num2str(session_num), ...
            ['scalp_freq_rtf_rc' zero_event '_Tafdf' num2str(session_num)],...
            ['condition_' conditions{cond_idx} '.nii']);
        V=spm_vol(img_path);
        for i=1:length(V)
            smooth_files{end+1,1}=[img_path ',' num2str(i)];
        end
        img_path=fullfile(params.data_dir, 'analysis', ...
            subj_info.subj_id, num2str(session_num), ...
            ['scalp_time_rtf_rc' zero_event '_Tafdf' num2str(session_num)],...
            ['condition_' conditions{cond_idx} '.nii']);
        V=spm_vol(img_path);
        for i=1:length(V)
            smooth_files{end+1,1}=[img_path ',' num2str(i)];
        end
        img_path=fullfile(params.data_dir, 'analysis', ...
            subj_info.subj_id, num2str(session_num), ...
            ['time_freq_rtf_rc' zero_event '_Tafdf' num2str(session_num)],...
            ['condition_' conditions{cond_idx} '.nii']);
        V=spm_vol(img_path);
        for i=1:length(V)
            smooth_files{end+1,1}=[img_path ',' num2str(i)];
        end
    end
    jobfile = 'smooth_job.m';
    inputs = {};
    inputs{1,1}=smooth_files;
    spm('defaults', 'EEG');
    spm_jobman('run', jobfile, inputs{:});
end

if params.run_stats
    jobfile = 'factorial_job.m';
    
    % Scalp x freq
    inputs={};
    inputs{1,1}={fullfile(params.data_dir, 'analysis', ...
            subj_info.subj_id, num2str(session_num), ...
            ['scalp_freq_rtf_rc' zero_event '_Tafdf' num2str(session_num)])};    
    for cond_idx=1:length(conditions)
        level_scans={};
        img_path=fullfile(params.data_dir, 'analysis', ...
            subj_info.subj_id, num2str(session_num), ...
            ['scalp_freq_rtf_rc' zero_event '_Tafdf' num2str(session_num)],...
            ['scondition_' conditions{cond_idx} '.nii']);
        V=spm_vol(img_path);
        for i=1:length(V)
            level_scans{end+1,1}=[img_path ',' num2str(i)];
        end
        inputs{1+cond_idx,1}=level_scans;
    end
    inputs{8,1}=3;
    inputs{9,1}=3;
    inputs{10,1}=3;
    spm('defaults', 'EEG');
    spm_jobman('run', jobfile, inputs{:});
    
    % Scalp x time
    inputs={};
    inputs{1,1}={fullfile(params.data_dir, 'analysis', ...
            subj_info.subj_id, num2str(session_num), ...
            ['scalp_time_rtf_rc' zero_event '_Tafdf' num2str(session_num)])};    
    for cond_idx=1:length(conditions)
        level_scans={};
        img_path=fullfile(params.data_dir, 'analysis', ...
            subj_info.subj_id, num2str(session_num), ...
            ['scalp_time_rtf_rc' zero_event '_Tafdf' num2str(session_num)],...
            ['scondition_' conditions{cond_idx} '.nii']);
        V=spm_vol(img_path);
        for i=1:length(V)
            level_scans{end+1,1}=[img_path ',' num2str(i)];
        end
        inputs{1+cond_idx,1}=level_scans;
    end
    inputs{8,1}=2;
    inputs{9,1}=2;
    inputs{10,1}=2;
    spm('defaults', 'EEG');
    spm_jobman('run', jobfile, inputs{:});
    
    % Time x frequency
    inputs={};
    inputs{1,1}={fullfile(params.data_dir, 'analysis', ...
            subj_info.subj_id, num2str(session_num), ...
            ['time_freq_rtf_rc' zero_event '_Tafdf' num2str(session_num)])};    
    for cond_idx=1:length(conditions)
        level_scans={};
        img_path=fullfile(params.data_dir, 'analysis', ...
            subj_info.subj_id, num2str(session_num), ...
            ['time_freq_rtf_rc' zero_event '_Tafdf' num2str(session_num)],...
            ['scondition_' conditions{cond_idx} '.nii']);
        V=spm_vol(img_path);
        for i=1:length(V)
            level_scans{end+1,1}=[img_path ',' num2str(i)];
        end
        inputs{1+cond_idx,1}=level_scans;
    end
    inputs{8,1}=4;
    inputs{9,1}=4;
    inputs{10,1}=4;
    spm('defaults', 'EEG');
    spm_jobman('run', jobfile, inputs{:});
end
