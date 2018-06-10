function plot_tfs(subj_info, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

subj_dir=fullfile('C:\meg_laminar\derivatives\spm12\', subj_info.subj_id);

conditions={'Undefined'};

dots_run_mean_tfs=[];
instr_run_mean_tfs=[];
resp_run_mean_tfs=[];
session_ids=[];

dots_session_mean_tfs=[];
instr_session_mean_tfs=[];
resp_session_mean_tfs=[];

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
    dots_session_mean_tfs(:,:,session_num)=mean(session_tfs,3);
    dots_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,1:round(size(session_tfs,3)/3)),3);
    run_idx=run_idx+1;
    dots_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,round(size(session_tfs,3)/3)+1:2*round(size(session_tfs,3)/3)),3);
    run_idx=run_idx+1;
    dots_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,2*round(size(session_tfs,3)/3)+1:end),3);
    run_idx=run_idx+1;
    session_ids(end+1)=session_num;
    session_ids(end+1)=session_num;
    session_ids(end+1)=session_num;
end

run_ICCs=[];
for session_num=1:length(subj_info.sessions)
    session_tfs=dots_run_mean_tfs(:,:,(session_num-1)*3+1:session_num*3);
    session_tfs=reshape(session_tfs,size(session_tfs,1)*size(session_tfs,2),size(session_tfs,3));
    run_ICCs(session_num)=IPN_icc(session_tfs,2,'k');    
end
session_ICC = IPN_icc(reshape(dots_session_mean_tfs,size(dots_session_mean_tfs,1)*size(dots_session_mean_tfs,2),size(dots_session_mean_tfs,3)),2,'k');
disp(sprintf('Dots - within-session mean ICC=%.4f, between-session ICC=%.4f', mean(run_ICCs), session_ICC));


clims=[min(dots_session_mean_tfs(:)) max(dots_session_mean_tfs(:))];
fig=figure('Position',[1 1 1600 300],'PaperUnits','points',...
    'PaperPosition',[1 1 1600 300],'PaperPositionMode','manual');
for session_num=1:length(subj_info.sessions)
    ax=subplot(1,length(subj_info.sessions),session_num);
    imagesc(times, freqs, squeeze(dots_session_mean_tfs(:,:,session_num)));%, clims);
    set(gca,'ydir','normal');
    colorbar();
    xlabel('Time');
    ylabel('Freq');
end

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
    instr_session_mean_tfs(:,:,session_num)=mean(session_tfs,3);
    instr_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,1:round(size(session_tfs,3)/3)),3);
    run_idx=run_idx+1;
    instr_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,round(size(session_tfs,3)/3)+1:2*round(size(session_tfs,3)/3)),3);
    run_idx=run_idx+1;
    instr_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,2*round(size(session_tfs,3)/3)+1:end),3);
    run_idx=run_idx+1;
end

run_ICCs=[];
for session_num=1:length(subj_info.sessions)
    session_tfs=instr_run_mean_tfs(:,:,(session_num-1)*3+1:session_num*3);
    session_tfs=reshape(session_tfs,size(session_tfs,1)*size(session_tfs,2),size(session_tfs,3));
    run_ICCs(session_num)=IPN_icc(session_tfs,2,'k');    
end
session_ICC = IPN_icc(reshape(instr_session_mean_tfs,size(instr_session_mean_tfs,1)*size(instr_session_mean_tfs,2),size(instr_session_mean_tfs,3)),2,'k');
disp(sprintf('Instr - within-session mean ICC=%.4f, between-session ICC=%.4f', mean(run_ICCs), session_ICC));

clims=[min(instr_session_mean_tfs(:)) max(instr_session_mean_tfs(:))];
fig=figure('Position',[1 1 1600 300],'PaperUnits','points',...
    'PaperPosition',[1 1 1600 300],'PaperPositionMode','manual');
for session_num=1:length(subj_info.sessions)
    ax=subplot(1,length(subj_info.sessions),session_num);
    imagesc(times, freqs, squeeze(instr_session_mean_tfs(:,:,session_num)));%, clims);
    set(gca,'ydir','normal');
    colorbar();
    xlabel('Time');
    ylabel('Freq');
end

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
    resp_session_mean_tfs(:,:,session_num)=mean(session_tfs,3);
    resp_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,1:round(size(session_tfs,3)/3)),3);
    run_idx=run_idx+1;
    resp_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,round(size(session_tfs,3)/3)+1:2*round(size(session_tfs,3)/3)),3);
    run_idx=run_idx+1;
    resp_run_mean_tfs(:,:,run_idx)=mean(session_tfs(:,:,2*round(size(session_tfs,3)/3)+1:end),3);
    run_idx=run_idx+1;
end

run_ICCs=[];
for session_num=1:length(subj_info.sessions)
    session_tfs=resp_run_mean_tfs(:,:,(session_num-1)*3+1:session_num*3);
    session_tfs=reshape(session_tfs,size(session_tfs,1)*size(session_tfs,2),size(session_tfs,3));
    run_ICCs(session_num)=IPN_icc(session_tfs,2,'k');    
end
session_ICC = IPN_icc(reshape(resp_session_mean_tfs,size(resp_session_mean_tfs,1)*size(resp_session_mean_tfs,2),size(resp_session_mean_tfs,3)),2,'k');
disp(sprintf('Resp - within-session mean ICC=%.4f, between-session ICC=%.4f', mean(run_ICCs), session_ICC));

clims=[min(resp_session_mean_tfs(:)) max(resp_session_mean_tfs(:))];
fig=figure('Position',[1 1 1600 300],'PaperUnits','points',...
    'PaperPosition',[1 1 1600 300],'PaperPositionMode','manual');
for session_num=1:length(subj_info.sessions)
    ax=subplot(1,length(subj_info.sessions),session_num);
    imagesc(times, freqs, squeeze(resp_session_mean_tfs(:,:,session_num)));%, clims);
    set(gca,'ydir','normal');
    colorbar();
    xlabel('Time');
    ylabel('Freq');
end

