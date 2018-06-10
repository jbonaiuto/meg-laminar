function export_data(subjects)

subj_info=subjects(3);
subj_dir=fullfile('D:\meg_laminar\derivatives\spm12\', subj_info.subj_id);

fid=fopen(fullfile(subj_dir, 'reproduciblity_data_topo.csv','w'));
fprintf(fid, 'Session,Run,AlignEvent,Channel,FieldIntensity\n');

dots_t=-2350;
instr_t=150;
resp_t=35;

instr_run_scalp_maps=[];
dots_run_scalp_maps=[];
resp_run_scalp_maps=[];

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
        instr_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),instr_t_idx,trials(1):trials(round(length(trials)/3))),3));
        instr_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),instr_t_idx,trials(round(length(trials)/3)+1):trials(2*round(length(trials)/3))),3));
        instr_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),instr_t_idx,trials(2*round(length(trials)/3)+1):trials(end)),3));
    else
        % Position of each meg channel
        ch_pos=instr_data.coor2D(meg_ch_idx);
        % Label for each meg channel
        ch_labels=instr_data.chanlabels(meg_ch_idx);
        instr_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,instr_t_idx,trials(1):trials(round(length(trials)/3))),3));
        instr_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,instr_t_idx,trials(round(length(trials)/3)+1):trials(2*round(length(trials)/3))),3));
        instr_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,instr_t_idx,trials(2*round(length(trials)/3)+1):trials(end)),3));
    end
    
    if session_num==1
        dots_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),dots_t_idx,trials(1):trials(round(length(trials)/3))),3));
        dots_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),dots_t_idx,trials(round(length(trials)/3)+1):trials(2*round(length(trials)/3))),3));
        dots_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx([1:70 72:275]),dots_t_idx,trials(2*round(length(trials)/3)+1):trials(end)),3));
    else
        dots_mean_scalp_vals=squeeze(mean(instr_data(meg_ch_idx,dots_t_idx,trials),3));
        dots_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,dots_t_idx,trials(1):trials(round(length(trials)/3))),3));
        dots_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,dots_t_idx,trials(round(length(trials)/3)+1):trials(2*round(length(trials)/3))),3));
        dots_run_scalp_maps(:,end+1)=squeeze(mean(instr_data(meg_ch_idx,dots_t_idx,trials(2*round(length(trials)/3)+1):trials(end)),3));
    end
        
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
        resp_run_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx([1:70 72:275]),resp_t_idx,trials(1):trials(round(length(trials)/3))),3));
        resp_run_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx([1:70 72:275]),resp_t_idx,trials(round(length(trials)/3)+1):trials(2*round(length(trials)/3))),3));
        resp_run_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx([1:70 72:275]),resp_t_idx,trials(2*round(length(trials)/3)+1):trials(end)),3));
    else
        % Position of each meg channel
        ch_pos=resp_data.coor2D(meg_ch_idx);
        % Label for each meg channel
        ch_labels=resp_data.chanlabels(meg_ch_idx);
        resp_run_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx,resp_t_idx,trials(1):trials(round(length(trials)/3))),3));
        resp_run_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx,resp_t_idx,trials(round(length(trials)/3)+1):trials(2*round(length(trials)/3))),3));
        resp_run_scalp_maps(:,end+1)=squeeze(mean(resp_data(meg_ch_idx,resp_t_idx,trials(2*round(length(trials)/3)+1):trials(end)),3));
    end    
end

r_idx=1;
for session_num=1:4
    for run=1:3
        for c_idx=1:length(ch_labels)
            fprintf(fid, sprintf('%d,%d,RDK,%s,%.4f\n', session_num, run, ch_labels{c_idx}, dots_run_scalp_maps(c_idx,r_idx)));
        end
        r_idx=r_idx+1;
    end
end

r_idx=1;
for session_num=1:4
    for run=1:3
        for c_idx=1:length(ch_labels)
            fprintf(fid, sprintf('%d,%d,InstructionCue,%s,%.4f\n', session_num, run, ch_labels{c_idx}, instr_run_scalp_maps(c_idx,r_idx)));
        end
        r_idx=r_idx+1;
    end
end

r_idx=1;
for session_num=1:4
    for run=1:3
        for c_idx=1:length(ch_labels)
            fprintf(fid, sprintf('%d,%d,Response,%s,%.4f\n', session_num, run, ch_labels{c_idx}, resp_run_scalp_maps(c_idx,r_idx)));
        end
        r_idx=r_idx+1;
    end
end

fclose(fid);


fid=fopen(fullfile(subj_dir, 'reproduciblity_data_erp.csv','w'));
fprintf(fid, 'Session,Run,AlignEvent,Time,FieldIntensity\n');

dots_ch='MLO31';
instr_ch='MLO32';
resp_ch='MLC17';
dots_min_time=-2.75;
dots_max_time=-1.75;
instr_min_time=-.25;
instr_max_time=.75;
resp_min_time=-.5;
resp_max_time=.5;

run_erps=[];
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
end

r_idx=1;
for session_num=1:4
    for run=1:3
        for t_idx=1:length(dots_time_idx)
            fprintf(fid,'%d,%d,RDK,%.4f,%.4f\n',session_num,run,time(dots_time_idx(t_idx)),run_erps(t_idx,r_idx));
        end
        r_idx=r_idx+1;
    end
end

run_erps=[];
for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, num2str(session_num));
    
    instr_file=fullfile(session_dir, sprintf('rcinstr_Tafdf%d.mat', session_num));
    instr_data=spm_eeg_load(instr_file);
    time=instr_data.time;
    
    instr_ch_idx=instr_data.indchannel(instr_ch);
    instr_trial_erps=squeeze(instr_data(instr_ch_idx,:,:));
    instr_time_idx=intersect(find(time>=instr_min_time),find(time<instr_max_time));
    
    run_erps(:,end+1)=mean(instr_trial_erps(instr_time_idx,1:size(instr_trial_erps,2)/3),2);
    run_erps(:,end+1)=mean(instr_trial_erps(instr_time_idx,size(instr_trial_erps,2)/3+1:2*size(instr_trial_erps,2)/3),2);
    run_erps(:,end+1)=mean(instr_trial_erps(instr_time_idx,2*size(instr_trial_erps,2)/3+1:size(instr_trial_erps,2)),2);
end

r_idx=1;
for session_num=1:4
    for run=1:3
        for t_idx=1:length(instr_time_idx)
            fprintf(fid,'%d,%d,InstructionCue,%.4f,%.4f\n',session_num,run,time(instr_time_idx(t_idx)),run_erps(t_idx,r_idx));
        end
        r_idx=r_idx+1;
    end
end

run_erps=[];
for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, num2str(session_num));
    
    resp_file=fullfile(session_dir, sprintf('rcresp_Tafdf%d.mat', session_num));
    resp_data=spm_eeg_load(resp_file);
    time=resp_data.time;
    
    resp_ch_idx=resp_data.indchannel(resp_ch);
    resp_trial_erps=squeeze(resp_data(resp_ch_idx,:,:));
    resp_time_idx=intersect(find(time>=resp_min_time),find(time<resp_max_time));
    
    run_erps(:,end+1)=mean(resp_trial_erps(resp_time_idx,1:size(resp_trial_erps,2)/3),2);
    run_erps(:,end+1)=mean(resp_trial_erps(resp_time_idx,size(resp_trial_erps,2)/3+1:2*size(resp_trial_erps,2)/3),2);
    run_erps(:,end+1)=mean(resp_trial_erps(resp_time_idx,2*size(resp_trial_erps,2)/3+1:size(resp_trial_erps,2)),2);
end

r_idx=1;
for session_num=1:4
    for run=1:3
        for t_idx=1:length(resp_time_idx)
            fprintf(fid,'%d,%d,Response,%.4f,%.4f\n',session_num,run,time(resp_time_idx(t_idx)),run_erps(t_idx,r_idx));
        end
        r_idx=r_idx+1;
    end
end

fclose(fid);


fid=fopen(fullfile(subj_dir, 'reproduciblity_data_tf_rdk.csv','w'));
fprintf(fid, 'Session,Run,AlignEvent,Time,Frequency,Power\n');

conditions={'Undefined'};

dots_run_mean_tfs=[];
instr_run_mean_tfs=[];
resp_run_mean_tfs=[];

run_idx=1;
for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-%02d',session_num));
    tf_dir=fullfile(session_dir,sprintf('rdots_tf_ffrcinstr_Tafdf%d',session_num));
    session_tfs=[];
    for cond_idx=1:length(conditions)
        condition=conditions{cond_idx};
        
        X=spm_vol(fullfile(tf_dir, sprintf('scondition_%s.nii', condition)));
        max_dim=max(X(1).dim);
        coords=X(1).mat*[[1:max_dim]' [1:max_dim]' ones(max_dim,1) ones(max_dim,1)]';
        
        times=coords(2,1:X(1).dim(2));
        freqs=coords(1,1:X(1).dim(1));
        
        img=spm_read_vols(X);
        img=squeeze(img(:,:,1,:));
        session_tfs(:,:,end+1:end+size(img,3))=img;
    end
    dots_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,1:round(size(session_tfs,3)/3)),3);
    run_idx=run_idx+1;
    dots_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,round(size(session_tfs,3)/3)+1:2*round(size(session_tfs,3)/3)),3);
    run_idx=run_idx+1;
    dots_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,2*round(size(session_tfs,3)/3)+1:end),3);
    run_idx=run_idx+1;
end

r_idx=1;
for session_num=1:4
    for run=1:3
        for t_idx=1:length(times)
            for f_idx=1:length(freqs)
                fprintf(fid,'%d,%d,RDK,%.4f,%.4f,%.4f\n', session_num,run,times(t_idx),freqs(f_idx),dots_run_mean_tfs(f_idx,t_idx,r_idx));
            end
        end
        r_idx=r_idx+1;
    end
end
fclose(fid);

fid=fopen(fullfile(subj_dir, 'reproduciblity_data_tf_instructioncue.csv','w'));
fprintf(fid, 'Session,Run,AlignEvent,Time,Frequency,Power\n');

run_idx=1;
for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-%02d',session_num));
    tf_dir=fullfile(session_dir,sprintf('rinstr_tf_ffrcinstr_Tafdf%d',session_num));
    session_tfs=[];
    for cond_idx=1:length(conditions)
        condition=conditions{cond_idx};
        
        X=spm_vol(fullfile(tf_dir, sprintf('scondition_%s.nii', condition)));
        max_dim=max(X(1).dim);
        coords=X(1).mat*[[1:max_dim]' [1:max_dim]' ones(max_dim,1) ones(max_dim,1)]';
        
        times=coords(2,1:X(1).dim(2));
        freqs=coords(1,1:X(1).dim(1));
        
        img=spm_read_vols(X);
        img=squeeze(img(:,:,1,:));
        session_tfs(:,:,end+1:end+size(img,3))=img;
    end
    instr_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,1:round(size(session_tfs,3)/3)),3);
    run_idx=run_idx+1;
    instr_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,round(size(session_tfs,3)/3)+1:2*round(size(session_tfs,3)/3)),3);
    run_idx=run_idx+1;
    instr_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,2*round(size(session_tfs,3)/3)+1:end),3);
    run_idx=run_idx+1;
end

r_idx=1;
for session_num=1:4
    for run=1:3
        for t_idx=1:length(times)
            for f_idx=1:length(freqs)
                fprintf(fid,'%d,%d,InstructionCue,%.4f,%.4f,%.4f\n', session_num,run,times(t_idx),freqs(f_idx),instr_run_mean_tfs(f_idx,t_idx,r_idx));
            end
        end
        r_idx=r_idx+1;
    end
end
fclose(fid);

fid=fopen(fullfile(subj_dir, 'reproduciblity_data_tf_resp.csv','w'));
fprintf(fid, 'Session,Run,AlignEvent,Time,Frequency,Power\n');

run_idx=1;
for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-%02d',session_num));
    tf_dir=fullfile(session_dir,sprintf('rresp_tf_ffrcresp_Tafdf%d',session_num));
    session_tfs=[];
    for cond_idx=1:length(conditions)
        condition=conditions{cond_idx};
        
        X=spm_vol(fullfile(tf_dir, sprintf('scondition_%s.nii', condition)));
        max_dim=max(X(1).dim);
        coords=X(1).mat*[[1:max_dim]' [1:max_dim]' ones(max_dim,1) ones(max_dim,1)]';
        
        times=coords(2,1:X(1).dim(2));
        freqs=coords(1,1:X(1).dim(1));
        
        img=spm_read_vols(X);
        img=squeeze(img(:,:,1,:));
        session_tfs(:,:,end+1:end+size(img,3))=img;
    end
    resp_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,1:round(size(session_tfs,3)/3)),3);
    run_idx=run_idx+1;
    resp_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,round(size(session_tfs,3)/3)+1:2*round(size(session_tfs,3)/3)),3);
    run_idx=run_idx+1;
    resp_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,2*round(size(session_tfs,3)/3)+1:end),3);
    run_idx=run_idx+1;
end

r_idx=1;
for session_num=1:4
    for run=1:3
        for t_idx=1:length(times)
            for f_idx=1:length(freqs)
                fprintf(fid,'%d,%d,Response,%.4f,%.4f,%.4f\n', session_num,run,times(t_idx),freqs(f_idx),resp_run_mean_tfs(f_idx,t_idx,r_idx));
            end
        end
        r_idx=r_idx+1;
    end
end
fclose(fid);