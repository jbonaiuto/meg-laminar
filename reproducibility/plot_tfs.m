function plot_tfs(subj_info, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

conditions={'Undefined'};

dots_run_mean_tfs=[];
instr_run_mean_tfs=[];
resp_run_mean_tfs=[];
session_ids=[];
subj_dir=fullfile('D:\pred_coding\derivatives\spm12', subj_info.subj_id);


dots_session_mean_tfs=[];
instr_session_mean_tfs=[];
resp_session_mean_tfs=[];

run_idx=1;
for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-0%d',session_num));
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
    session_dir=fullfile(subj_dir, sprintf('ses-0%d',session_num));
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
    session_dir=fullfile(subj_dir, sprintf('ses-0%d',session_num));
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


run_all_r_sq=[];
within_r_sq=[];
for i=1:size(instr_run_mean_tfs,3)
    map_one=squeeze(instr_run_mean_tfs(:,:,i));
    for j=1:size(instr_run_mean_tfs,3)
        map_two=squeeze(instr_run_mean_tfs(:,:,j));
        mdl=fitlm(map_one(:),map_two(:));
        run_all_r_sq(i,j)=mdl.Rsquared.Adjusted;
        if i<j
            if session_ids(i)==session_ids(j)
                within_r_sq(end+1)=mdl.Rsquared.Adjusted;
            end
        end
    end
end
session_all_r_sq=[];
between_r_sq=[];
for i=1:size(instr_session_mean_tfs,3)
    map_one=squeeze(instr_session_mean_tfs(:,:,i));
    for j=1:size(instr_session_mean_tfs,3)
        map_two=squeeze(instr_session_mean_tfs(:,:,j));
        mdl=fitlm(map_one(:),map_two(:));
        session_all_r_sq(i,j)=mdl.Rsquared.Adjusted;
        if i<j
            between_r_sq(end+1)=mdl.Rsquared.Adjusted;
        end
    end
end

figure();
subplot(1,2,1);
imagesc(run_all_r_sq);
set(gca,'clim',[0 1]);
colorbar();
axis square;
set(gca,'xtick',[2 5 8 11]);
set(gca,'xticklabel',{'Session 1','Session 2','Session 3','Session 4'});
xticklabel_rotate([],45);
set(gca,'ytick',[2 5 8 11]);
set(gca,'yticklabel',{'Session 1','Session 2','Session 3','Session 4'});

subplot(1,2,2);
imagesc(session_all_r_sq);
set(gca,'clim',[0 1]);
colorbar();
axis square;
set(gca,'xtick',[1 2 3 4]);
set(gca,'xticklabel',{'Session 1','Session 2','Session 3','Session 4'});
xticklabel_rotate([],45);
set(gca,'ytick',[1 2 3 4]);
set(gca,'yticklabel',{'Session 1','Session 2','Session 3','Session 4'});

disp(sprintf('Instr - within-session mean R^2=%.4f, between-session mean R^2=%.4f', mean(within_r_sq), mean(between_r_sq)));




run_all_r_sq=[];
within_r_sq=[];
for i=1:size(dots_run_mean_tfs,3)
    map_one=squeeze(dots_run_mean_tfs(:,:,i));
    for j=1:size(dots_run_mean_tfs,3)
        map_two=squeeze(dots_run_mean_tfs(:,:,j));
        mdl=fitlm(map_one(:),map_two(:));
        run_all_r_sq(i,j)=mdl.Rsquared.Adjusted;
        if i<j
            if session_ids(i)==session_ids(j)
                within_r_sq(end+1)=mdl.Rsquared.Adjusted;
            end
        end
    end
end
session_all_r_sq=[];
between_r_sq=[];
for i=1:size(dots_session_mean_tfs,3)
    map_one=squeeze(dots_session_mean_tfs(:,:,i));
    for j=1:size(dots_session_mean_tfs,3)
        map_two=squeeze(dots_session_mean_tfs(:,:,j));
        mdl=fitlm(map_one(:),map_two(:));
        session_all_r_sq(i,j)=mdl.Rsquared.Adjusted;
        if i<j
            between_r_sq(end+1)=mdl.Rsquared.Adjusted;
        end
    end
end

figure();
subplot(1,2,1);
imagesc(run_all_r_sq);
set(gca,'clim',[0 1]);
colorbar();
axis square;
set(gca,'xtick',[2 5 8 11]);
set(gca,'xticklabel',{'Session 1','Session 2','Session 3','Session 4'});
xticklabel_rotate([],45);
set(gca,'ytick',[2 5 8 11]);
set(gca,'yticklabel',{'Session 1','Session 2','Session 3','Session 4'});

subplot(1,2,2);
imagesc(session_all_r_sq);
set(gca,'clim',[0 1]);
colorbar();
axis square;
set(gca,'xtick',[1 2 3 4]);
set(gca,'xticklabel',{'Session 1','Session 2','Session 3','Session 4'});
xticklabel_rotate([],45);
set(gca,'ytick',[1 2 3 4]);
set(gca,'yticklabel',{'Session 1','Session 2','Session 3','Session 4'});

disp(sprintf('Dots - within-session mean R^2=%.4f, between-session mean R^2=%.4f', mean(within_r_sq), mean(between_r_sq)));



run_all_r_sq=[];
within_r_sq=[];
for i=1:size(resp_run_mean_tfs,3)
    map_one=squeeze(resp_run_mean_tfs(:,:,i));
    for j=1:size(resp_run_mean_tfs,3)
        map_two=squeeze(resp_run_mean_tfs(:,:,j));
        mdl=fitlm(map_one(:),map_two(:));
        run_all_r_sq(i,j)=mdl.Rsquared.Adjusted;
        if i<j
            if session_ids(i)==session_ids(j)
                within_r_sq(end+1)=mdl.Rsquared.Adjusted;
            end
        end
    end
end
session_all_r_sq=[];
between_r_sq=[];
for i=1:size(resp_session_mean_tfs,3)
    map_one=squeeze(resp_session_mean_tfs(:,:,i));
    for j=1:size(resp_session_mean_tfs,3)
        map_two=squeeze(resp_session_mean_tfs(:,:,j));
        mdl=fitlm(map_one(:),map_two(:));
        session_all_r_sq(i,j)=mdl.Rsquared.Adjusted;
        if i<j
            between_r_sq(end+1)=mdl.Rsquared.Adjusted;
        end
    end
end

figure();
subplot(1,2,1);
imagesc(run_all_r_sq);
set(gca,'clim',[0 1]);
colorbar();
axis square;
set(gca,'xtick',[2 5 8 11]);
set(gca,'xticklabel',{'Session 1','Session 2','Session 3','Session 4'});
xticklabel_rotate([],45);
set(gca,'ytick',[2 5 8 11]);
set(gca,'yticklabel',{'Session 1','Session 2','Session 3','Session 4'});

subplot(1,2,2);
imagesc(session_all_r_sq);
set(gca,'clim',[0 1]);
colorbar();
axis square;
set(gca,'xtick',[1 2 3 4]);
set(gca,'xticklabel',{'Session 1','Session 2','Session 3','Session 4'});
xticklabel_rotate([],45);
set(gca,'ytick',[1 2 3 4]);
set(gca,'yticklabel',{'Session 1','Session 2','Session 3','Session 4'});

disp(sprintf('Resp - within-session mean R^2=%.4f, between-session mean R^2=%.4f', mean(within_r_sq), mean(between_r_sq)));
