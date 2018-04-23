function plot_field_maps(subj_info, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

subj_dir=fullfile('D:\pred_coding\analysis\', subj_info.subj_id);
%subj_dir=fullfile('D:\pred_coding\derivatives\spm12', subj_info.subj_id);

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

for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, num2str(session_num));
%    session_dir=fullfile(subj_dir, sprintf('ses-0%d',session_num));
    
    instr_file=fullfile(session_dir, sprintf('rcinstr_Tafdf%d.mat', session_num));
    instr_data=spm_eeg_load(instr_file);
    % All meg channel idx
    meg_ch_idx = instr_data.indchantype('MEG');        
    % Position of each meg channel
    ch_pos=instr_data.coor2D(meg_ch_idx);
    % Label for each meg channel
    ch_labels=instr_data.chanlabels(meg_ch_idx);

    instr_t_idx=min(find(instr_data.time([],'ms')>=instr_t));
    dots_t_idx=min(find(instr_data.time([],'ms')>=dots_t));
    trials=setdiff([1:size(instr_data,3)],instr_data.badtrials);
    
    instr_mean_scalp_vals=squeeze(mean(instr_data(meg_ch_idx,instr_t_idx,trials),3));
    dots_mean_scalp_vals=squeeze(mean(instr_data(meg_ch_idx,dots_t_idx,trials),3));
    
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
    % Position of each meg channel
    ch_pos=resp_data.coor2D(meg_ch_idx);
    % Label for each meg channel
    ch_labels=resp_data.chanlabels(meg_ch_idx);

    resp_t_idx=min(find(resp_data.time([],'ms')>=resp_t));
    trials=setdiff([1:size(resp_data,3)],resp_data.badtrials);
    
    resp_mean_scalp_vals=squeeze(mean(resp_data(meg_ch_idx,resp_t_idx,trials),3));
    
    sessions_resp(session_num).scalp_vals=resp_mean_scalp_vals;
    sessions_resp(session_num).ch_pos=ch_pos;
    sessions_resp(session_num).ch_labels=ch_labels;
    
    resp_min_val=min([resp_min_val, min(resp_mean_scalp_vals(:))]);
    resp_max_val=max([resp_max_val, max(resp_mean_scalp_vals(:))]);
end

fig=figure('Position',[1 1 1600 400],'PaperUnits','points',...
    'PaperPosition',[1 1 1600 400],'PaperPositionMode','manual');
for session_num=1:length(subj_info.sessions)
    ax=subplot(1,length(subj_info.sessions),session_num);
    in.f=fig;
    in.ParentAxes=ax;
    in.noButtons=true;
    in.type='MEG';
    in.min=instr_min_val;
    in.max=instr_max_val;
    [ZI,f]=spm_eeg_plotScalpData(sessions_instr(session_num).scalp_vals,sessions_instr(session_num).ch_pos,sessions_instr(session_num).ch_labels,in);
end

fig=figure('Position',[1 1 1600 400],'PaperUnits','points','PaperPosition',[1 1 1600 400],'PaperPositionMode','manual');
for session_num=1:length(subj_info.sessions)
    ax=subplot(1,length(subj_info.sessions),session_num);
    in.f=fig;
    in.ParentAxes=ax;
    in.noButtons=true;
    in.type='MEG';
    in.min=dots_min_val;
    in.max=dots_max_val;
    [ZI,f]=spm_eeg_plotScalpData(sessions_dots(session_num).scalp_vals,sessions_dots(session_num).ch_pos,sessions_dots(session_num).ch_labels,in);
end

fig=figure('Position',[1 1 1600 400],'PaperUnits','points','PaperPosition',[1 1 1600 400],'PaperPositionMode','manual');
for session_num=1:length(subj_info.sessions)
    ax=subplot(1,length(subj_info.sessions),session_num);
    in.f=fig;
    in.ParentAxes=ax;
    in.noButtons=true;
    in.type='MEG';
    in.min=resp_min_val;
    in.max=resp_max_val;
    [ZI,f]=spm_eeg_plotScalpData(sessions_resp(session_num).scalp_vals,sessions_resp(session_num).ch_pos,sessions_resp(session_num).ch_labels,in);
end
