function plot_tfs(subj_info, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

subj_dir=fullfile('C:\pred_coding\analysis\', subj_info.subj_id);

conditions={'congruent-low','congruent-med','congruent-high',...
    'incongruent-low','incongruent-med','incongruent-high'};

dots_session_mean_tfs=[];
instr_session_mean_tfs=[];
resp_session_mean_tfs=[];

for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, num2str(session_num));
    tf_dir=fullfile(session_dir,sprintf('rtf_dotstf_rcinstr_Tafdf%d',session_num));
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
end

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

for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, num2str(session_num));
    tf_dir=fullfile(session_dir,sprintf('rtf_instrtf_rcinstr_Tafdf%d',session_num));
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
end

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

for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, num2str(session_num));
    tf_dir=fullfile(session_dir,sprintf('rtf_rcresp_Tafdf%d',session_num));
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
end

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