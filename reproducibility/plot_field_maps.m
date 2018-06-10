function plot_field_maps(subj_info, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

subj_dir=fullfile('D:\meg_laminar\derivatives\spm12', subj_info.subj_id);

dots_t=-2350;
instr_t=150;
resp_t=35;

sessions_instr=[];
sessions_dots=[];
sessions_resp=[];

instr_min_val=inf;
instr_max_val=-inf;
dots_min_val=inf;
dots_max_val=-inf;
resp_min_val=inf;
resp_max_val=-inf;

instr_run_scalp_maps=[];
dots_run_scalp_maps=[];
resp_run_scalp_maps=[];

instr_session_scalp_maps=[];
dots_session_scalp_maps=[];
resp_session_scalp_maps=[];

session_ids=[];

for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-0%d',session_num));
    
    instr_file=fullfile(session_dir, sprintf('rcinstr_Tafdf%d.mat', session_num));
    instr_data=spm_eeg_load(instr_file);
    % All meg channel idx
    meg_ch_idx = instr_data.indchantype('MEG');        
    

    instr_t_idx=min(find(instr_data.time([],'ms')>=instr_t));
    dots_t_idx=min(find(instr_data.time([],'ms')>=dots_t));
    trials=setdiff([1:size(instr_data,3)],instr_data.badtrials);
    
    
    if session_num==1
        % Position of each meg channel
        ch_pos=instr_data.coor2D(meg_ch_idx([1:70 72:275]));
        % Label for each meg channel
        ch_labels=instr_data.chanlabels(meg_ch_idx([1:70 72:275]));
        instr_mean_scalp_vals=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),instr_t_idx,trials),3));
        instr_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),instr_t_idx,trials(1):trials(round(length(trials)/3))),3));
        instr_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),instr_t_idx,trials(round(length(trials)/3)+1):trials(2*round(length(trials)/3))),3));
        instr_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),instr_t_idx,trials(2*round(length(trials)/3)+1):trials(end)),3));
        instr_session_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),instr_t_idx,trials),3));
    else
        % Position of each meg channel
        ch_pos=instr_data.coor2D(meg_ch_idx);
        % Label for each meg channel
        ch_labels=instr_data.chanlabels(meg_ch_idx);
        instr_mean_scalp_vals=squeeze(mean(instr_data(meg_ch_idx,instr_t_idx,trials),3));
        instr_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,instr_t_idx,trials(1):trials(round(length(trials)/3))),3));
        instr_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,instr_t_idx,trials(round(length(trials)/3)+1):trials(2*round(length(trials)/3))),3));
        instr_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,instr_t_idx,trials(2*round(length(trials)/3)+1):trials(end)),3));
        instr_session_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,instr_t_idx,trials),3));
    end
    
    if session_num==1
        dots_mean_scalp_vals=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),dots_t_idx,trials),3));
        dots_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),dots_t_idx,trials(1):trials(round(length(trials)/3))),3));
        dots_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),dots_t_idx,trials(round(length(trials)/3)+1):trials(2*round(length(trials)/3))),3));
        dots_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),dots_t_idx,trials(2*round(length(trials)/3)+1):trials(end)),3));
        dots_session_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),dots_t_idx,trials),3));
    else
        dots_mean_scalp_vals=squeeze(mean(instr_data(meg_ch_idx,dots_t_idx,trials),3));
        dots_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,dots_t_idx,trials(1):trials(round(length(trials)/3))),3));
        dots_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,dots_t_idx,trials(round(length(trials)/3)+1):trials(2*round(length(trials)/3))),3));
        dots_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,dots_t_idx,trials(2*round(length(trials)/3)+1):trials(end)),3));
        dots_session_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,dots_t_idx,trials),3));
    end
    
    sessions_instr(session_num).scalp_vals=instr_mean_scalp_vals;
    sessions_instr(session_num).ch_pos=ch_pos;
    sessions_instr(session_num).ch_labels=ch_labels;
    
    sessions_dots(session_num).scalp_vals=dots_mean_scalp_vals;
    sessions_dots(session_num).ch_pos=ch_pos;
    sessions_dots(session_num).ch_labels=ch_labels;
    
    instr_min_val=min([instr_min_val, min(instr_mean_scalp_vals(:))]);
    instr_max_val=max([instr_max_val, max(instr_mean_scalp_vals(:))]);
    
    dots_min_val=min([dots_min_val, min(dots_mean_scalp_vals(:))]);
    dots_max_val=max([dots_max_val, max(dots_mean_scalp_vals(:))]);
    
    resp_file=fullfile(session_dir, sprintf('rcresp_Tafdf%d.mat', session_num));
    resp_data=spm_eeg_load(resp_file);
    % All meg channel idx
    meg_ch_idx = resp_data.indchantype('MEG');        
    

    resp_t_idx=min(find(resp_data.time([],'ms')>=resp_t));
    trials=setdiff([1:size(resp_data,3)],resp_data.badtrials);
    
    
    if session_num==1
        % Position of each meg channel
        ch_pos=resp_data.coor2D(meg_ch_idx([1:70 72:275]));
        % Label for each meg channel
        ch_labels=resp_data.chanlabels(meg_ch_idx([1:70 72:275]));
        resp_mean_scalp_vals=squeeze(mean(resp_data(meg_ch_idx([1:70 72:275]),resp_t_idx,trials),3));
        resp_run_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx([1:70 72:275]),resp_t_idx,trials(1):trials(round(length(trials)/3))),3));
        resp_run_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx([1:70 72:275]),resp_t_idx,trials(round(length(trials)/3)+1):trials(2*round(length(trials)/3))),3));
        resp_run_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx([1:70 72:275]),resp_t_idx,trials(2*round(length(trials)/3)+1):trials(end)),3));
        resp_session_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx([1:70 72:275]),resp_t_idx,trials),3));
    else
        % Position of each meg channel
        ch_pos=resp_data.coor2D(meg_ch_idx);
        % Label for each meg channel
        ch_labels=resp_data.chanlabels(meg_ch_idx);
        resp_mean_scalp_vals=squeeze(mean(resp_data(meg_ch_idx,resp_t_idx,trials),3));
        resp_run_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx,resp_t_idx,trials(1):trials(round(length(trials)/3))),3));
        resp_run_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx,resp_t_idx,trials(round(length(trials)/3)+1):trials(2*round(length(trials)/3))),3));
        resp_run_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx,resp_t_idx,trials(2*round(length(trials)/3)+1):trials(end)),3));
        resp_session_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx,resp_t_idx,trials),3));
    end
    session_ids(end+1)=session_num;
    session_ids(end+1)=session_num;
    session_ids(end+1)=session_num;
    
    sessions_resp(session_num).scalp_vals=resp_mean_scalp_vals;
    sessions_resp(session_num).ch_pos=ch_pos;
    sessions_resp(session_num).ch_labels=ch_labels;
    
    resp_min_val=min([resp_min_val, min(resp_mean_scalp_vals(:))]);
    resp_max_val=max([resp_max_val, max(resp_mean_scalp_vals(:))]);

end

run_ICCs=[];
for session_num=1:length(subj_info.sessions)
    session_scalp_maps=dots_run_scalp_maps(:,(session_num-1)*3+1:session_num*3);
    run_ICCs(session_num)=IPN_icc(session_scalp_maps,2,'k');    
end
session_ICC = IPN_icc(dots_session_scalp_maps,2,'k');
disp(sprintf('Dots - within-session mean ICC=%.4f, between-session ICC=%.4f', mean(run_ICCs), session_ICC));

run_ICCs=[];
for session_num=1:length(subj_info.sessions)
    session_scalp_maps=instr_run_scalp_maps(:,(session_num-1)*3+1:session_num*3);
    run_ICCs(session_num)=IPN_icc(session_scalp_maps,2,'k');    
end
session_ICC = IPN_icc(instr_session_scalp_maps,2,'k');
disp(sprintf('Instr - within-session mean ICC=%.4f, between-session ICC=%.4f', mean(run_ICCs), session_ICC));

run_ICCs=[];
for session_num=1:length(subj_info.sessions)
    session_scalp_maps=resp_run_scalp_maps(:,(session_num-1)*3+1:session_num*3);
    run_ICCs(session_num)=IPN_icc(session_scalp_maps,2,'k');    
end
session_ICC = IPN_icc(resp_session_scalp_maps,2,'k');
disp(sprintf('Resp - within-session mean ICC=%.4f, between-session ICC=%.4f', mean(run_ICCs), session_ICC));
