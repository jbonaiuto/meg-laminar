function plot_erps(subj_info, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

subj_dir=fullfile('D:\pred_coding\derivatives\spm12\', subj_info.subj_id);

dots_ch='MLO31';
instr_ch='MLO32';
resp_ch='MLC17';
dots_min_time=-2.75;
dots_max_time=-1.75;
instr_min_time=-.25;
instr_max_time=.75;
resp_min_time=-.5;
resp_max_time=.5;


figure();
hold on;
for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-0%d',session_num));
    
    instr_file=fullfile(session_dir, sprintf('rcinstr_Tafdf%d.mat', session_num));
    instr_data=spm_eeg_load(instr_file);
    time=instr_data.time;
    
    dots_ch_idx=instr_data.indchannel(dots_ch);
    dots_trial_erps=squeeze(instr_data(dots_ch_idx,:,:));
    dots_time_idx=intersect(find(time>=dots_min_time),find(time<dots_max_time));
    
    mean_erp=mean(dots_trial_erps,2);
    stderr_erp=std(dots_trial_erps,0,2)/sqrt(size(dots_trial_erps,2));
    H=shadedErrorBar(time(dots_time_idx)+2.5,mean_erp(dots_time_idx),stderr_erp(dots_time_idx),'b');
end
xlim([dots_min_time dots_max_time]+2.5);
xlabel('Time (s)');
ylabel('fT');

figure();
hold on;
for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-0%d',session_num));
    
    instr_file=fullfile(session_dir, sprintf('rcinstr_Tafdf%d.mat', session_num));
    instr_data=spm_eeg_load(instr_file);
    time=instr_data.time;
    
    instr_ch_idx=instr_data.indchannel(instr_ch);
    instr_trial_erps=squeeze(instr_data(instr_ch_idx,:,:));
    instr_time_idx=intersect(find(time>=instr_min_time),find(time<instr_max_time));
    
    mean_erp=mean(instr_trial_erps,2);
    stderr_erp=std(instr_trial_erps,0,2)/sqrt(size(instr_trial_erps,2));
    H=shadedErrorBar(time(instr_time_idx),mean_erp(instr_time_idx),stderr_erp(instr_time_idx),'b');
end
xlim([instr_min_time instr_max_time]);
xlabel('Time (s)');
ylabel('fT');

figure();
hold on;
for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-0%d',session_num));
    
    resp_file=fullfile(session_dir, sprintf('rcresp_Tafdf%d.mat', session_num));
    resp_data=spm_eeg_load(resp_file);
    time=resp_data.time;
    
    resp_ch_idx=resp_data.indchannel(resp_ch);
    resp_trial_erps=squeeze(resp_data(resp_ch_idx,:,:));
    resp_time_idx=intersect(find(time>=resp_min_time),find(time<resp_max_time));
    
    mean_erp=mean(resp_trial_erps,2);
    stderr_erp=std(resp_trial_erps,0,2)/sqrt(size(resp_trial_erps,2));
    H=shadedErrorBar(time(resp_time_idx),mean_erp(resp_time_idx),stderr_erp(resp_time_idx),'b');
end
xlim([resp_min_time resp_max_time]);
xlabel('Time (s)');
ylabel('fT');
