function preprocess_session(subj_info, session_num, run_num, varargin)

defaults = struct('delete_last_resp',false,'highpass_freq',2.0,'downsample',250,'lowpass_freq',100.0,'dots_epoch',[-1000 1000],'instr_epoch',[-3500 1500],'resp_epoch',[-2000 2000], 'remove_blinks', false, 'blink_channel','MLT31','manual_reject',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Directory containing session data
session_dir=fullfile('/data','pred_coding','scanning', subj_info.subj_id, num2str(session_num));
% Directory to put results
analysis_dir=fullfile('/data','pred_coding','analysis',subj_info.subj_id,num2str(session_num));
if exist(analysis_dir,'dir')~=7
    mkdir(analysis_dir);
end
% Directory to put report
report_dir=fullfile(analysis_dir,'preprocessing');
if exist(report_dir,'dir')~=7
    mkdir(report_dir);
end

% Add path to report generation code
addpath /home/jbonaiuto/Documents/MATLAB/m2html/
tpl = template('/home/jbonaiuto/Projects/meg/pred_coding/src/matlab/analysis_v2/preprocessing/templates/','keep');
tpl = set(tpl,'file',{'run'},{'run_template.tpl'});
tpl = set(tpl,'var',{'PAGETITLE'},{[subj_info.subj_id subj_info.birth_date ': Session ' num2str(session_num)]});

spm('defaults','eeg');
spm_jobman('initcfg');

for run_num=1:subj_info.sessions(session_num)
    % Directory containing run data
    run_code=[subj_info.subj_id subj_info.birth_date '_JamesBonaiuto_' subj_info.scan_date{session_num} '_0' num2str(run_num)];
    run_dir=fullfile(session_dir, [run_code '.ds']);

    % File containing stimulus information for each trial
    stim_file=fullfile(session_dir, ['stim_' subj_info.subj_id '_' num2str(run_num) '.mat']);
    % File containing behavioral data from each trial
    data_file=fullfile(session_dir, ['data_' subj_info.subj_id '_' num2str(run_num) '.mat']);

    clear jobs
    batch_idx=1;
    matlabbatch={};

    
end

