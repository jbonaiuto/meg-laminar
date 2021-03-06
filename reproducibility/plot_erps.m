function plot_erps(subj_info, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

subj_dir=fullfile('D:\meg_laminar\derivatives\spm12\', subj_info.subj_id);

dots_ch='MLO31';
instr_ch='MLO32';
resp_ch='MLC17';
dots_min_time=-2.75;
dots_max_time=-1.75;
instr_min_time=-.25;
instr_max_time=.75;
resp_min_time=-.5;
resp_max_time=.5;

session_colors=[102,194,165; 252,141,98; 141,160,203; 231,138,195]./255.0;
all_colors=[166,206,227; 31,120,180; 178,223,138; 51,160,44; 251,154,153; 227,26,28; 253,191,111; 255,127,0; 202,178,214; 106,61,154; 255,255,153; 177,89,40]./255.0;

figure();
hold on;
run_erps=[];
session_ids=[];
for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-%02d',session_num));
    
    instr_file=fullfile(session_dir, sprintf('rcinstr_Tafdf%d.mat', session_num));
    instr_data=spm_eeg_load(instr_file);
    time=instr_data.time;
    
    dots_ch_idx=instr_data.indchannel(dots_ch);
    dots_trial_erps=squeeze(instr_data(dots_ch_idx,:,:));
    dots_time_idx=intersect(find(time>=dots_min_time),find(time<dots_max_time));
    
    run_erps(:,end+1)=mean(dots_trial_erps(dots_time_idx,1:size(dots_trial_erps,2)/3),2);
    run_erps(:,end+1)=mean(dots_trial_erps(dots_time_idx,size(dots_trial_erps,2)/3+1:2*size(dots_trial_erps,2)/3),2);
    run_erps(:,end+1)=mean(dots_trial_erps(dots_time_idx,2*size(dots_trial_erps,2)/3+1:size(dots_trial_erps,2)),2);
    session_ids(1,end+1)=session_num;
    session_ids(1,end+1)=session_num;
    session_ids(1,end+1)=session_num;
    
    mean_erp=mean(dots_trial_erps,2);
    stderr_erp=std(dots_trial_erps,0,2)/sqrt(size(dots_trial_erps,2));
    H=shadedErrorBar(time(dots_time_idx)+2.5,mean_erp(dots_time_idx),stderr_erp(dots_time_idx),'b');
end
xlim([dots_min_time dots_max_time]+2.5);
xlabel('Time (s)');
ylabel('fT');

mean_session_erps=[];
run_ICCs=[];
for session_num=1:length(subj_info.sessions)
    run_ICCs(session_num)=IPN_icc(run_erps(:,(session_num-1)*3+1:session_num*3),2,'k');
    mean_session_erps(:,end+1)=mean(run_erps(:,(session_num-1)*3+1:session_num*3),2);
end
session_ICC = IPN_icc(mean_session_erps,2,'k');
disp(sprintf('Dots - within-session mean ICC=%.4f, between-session ICC=%.4f', mean(run_ICCs), session_ICC));

figure();
hold on;
run_erps=[];
session_ids=[];
for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-%02d',session_num));
    
    instr_file=fullfile(session_dir, sprintf('rcinstr_Tafdf%d.mat', session_num));
    instr_data=spm_eeg_load(instr_file);
    time=instr_data.time;
    
    instr_ch_idx=instr_data.indchannel(instr_ch);
    instr_trial_erps=squeeze(instr_data(instr_ch_idx,:,:));
    instr_time_idx=intersect(find(time>=instr_min_time),find(time<instr_max_time));
    
    run_erps(:,end+1)=mean(instr_trial_erps(instr_time_idx,1:size(instr_trial_erps,2)/3),2);
    run_erps(:,end+1)=mean(instr_trial_erps(instr_time_idx,size(instr_trial_erps,2)/3+1:2*size(instr_trial_erps,2)/3),2);
    run_erps(:,end+1)=mean(instr_trial_erps(instr_time_idx,2*size(instr_trial_erps,2)/3+1:size(instr_trial_erps,2)),2);
    session_ids(1,end+1)=session_num;
    session_ids(1,end+1)=session_num;
    session_ids(1,end+1)=session_num;
    
    mean_erp=mean(instr_trial_erps,2);
    stderr_erp=std(instr_trial_erps,0,2)/sqrt(size(instr_trial_erps,2));
    H=shadedErrorBar(time(instr_time_idx),mean_erp(instr_time_idx),stderr_erp(instr_time_idx),'b');
end
xlim([instr_min_time instr_max_time]);
xlabel('Time (s)');
ylabel('fT');

mean_session_erps=[];
run_ICCs=[];
for session_num=1:length(subj_info.sessions)
    run_ICCs(session_num)=IPN_icc(run_erps(:,(session_num-1)*3+1:session_num*3),2,'k');
    mean_session_erps(:,end+1)=mean(run_erps(:,(session_num-1)*3+1:session_num*3),2);
end
session_ICC = IPN_icc(mean_session_erps,2,'k');
disp(sprintf('Instr - within-session mean ICC=%.4f, between-session ICC=%.4f', mean(run_ICCs), session_ICC));

figure();
hold on;
run_erps=[];
session_ids=[];
for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-%02d',session_num));
    
    resp_file=fullfile(session_dir, sprintf('rcresp_Tafdf%d.mat', session_num));
    resp_data=spm_eeg_load(resp_file);
    time=resp_data.time;
    
    resp_ch_idx=resp_data.indchannel(resp_ch);
    resp_trial_erps=squeeze(resp_data(resp_ch_idx,:,:));
    resp_time_idx=intersect(find(time>=resp_min_time),find(time<resp_max_time));
    
    run_erps(:,end+1)=mean(resp_trial_erps(resp_time_idx,1:size(resp_trial_erps,2)/3),2);
    run_erps(:,end+1)=mean(resp_trial_erps(resp_time_idx,size(resp_trial_erps,2)/3+1:2*size(resp_trial_erps,2)/3),2);
    run_erps(:,end+1)=mean(resp_trial_erps(resp_time_idx,2*size(resp_trial_erps,2)/3+1:size(resp_trial_erps,2)),2);
    session_ids(1,end+1)=session_num;
    session_ids(1,end+1)=session_num;
    session_ids(1,end+1)=session_num;
    
    mean_erp=mean(resp_trial_erps,2);
    stderr_erp=std(resp_trial_erps,0,2)/sqrt(size(resp_trial_erps,2));
    H=shadedErrorBar(time(resp_time_idx),mean_erp(resp_time_idx),stderr_erp(resp_time_idx),'b');
end
xlim([resp_min_time resp_max_time]);
xlabel('Time (s)');
ylabel('fT');

mean_session_erps=[];
run_ICCs=[];
for session_num=1:length(subj_info.sessions)
    run_ICCs(session_num)=IPN_icc(run_erps(:,(session_num-1)*3+1:session_num*3),2,'k');
    mean_session_erps(:,end+1)=mean(run_erps(:,(session_num-1)*3+1:session_num*3),2);
end
session_ICC = IPN_icc(mean_session_erps,2,'k');
disp(sprintf('Resp - within-session mean ICC=%.4f, between-session ICC=%.4f', mean(run_ICCs), session_ICC));

