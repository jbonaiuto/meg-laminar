function average_subject_tfs(subj_info)

epochs={'dots','instr','resp'};
align_evts={'instr','instr','resp'};

spm('defaults','eeg');

subj_dir=fullfile('C:\mag_laminar\derivatives\spm12\', subj_info.subj_id);

sessions=[1:length(subj_info.sessions)];
if strcmp(subj_info.subj_id,'nc')
    sessions=[3];
end

for epoch_idx=1:length(epochs)
    epoch=epochs{epoch_idx};
    align_evt=align_evts{epoch_idx};
    
    epoch_dir=fullfile(subj_dir,epoch);
    mkdir(epoch_dir);
    
    session_tfs=[];
    session_smoothed_tfs=[];
    for session_idx=1:length(sessions)
        session_num=sessions(session_idx);
        session_dir=fullfile(subj_dir, sprintf('ses-%02d',session_num));
        
        tf_dir=fullfile(session_dir,sprintf('r%s_tf_ffrc%s_Tafdf%d',epoch,align_evt,session_num));
        
        X=spm_vol(fullfile(tf_dir, 'condition_Undefined.nii'));
        max_dim=max(X(1).dim);
        coords=X(1).mat*[[1:max_dim]' [1:max_dim]' ones(max_dim,1) ones(max_dim,1)]';

        times=coords(2,1:X(1).dim(2));
        freqs=coords(1,1:X(1).dim(1));

        img=spm_read_vols(X);
        img=squeeze(img(:,:,1,:));
        session_tfs(:,:,end+1:end+size(img,3))=img;
        
        sX=spm_vol(fullfile(tf_dir, 'scondition_Undefined.nii'));
        img=spm_read_vols(sX);
        img=squeeze(img(:,:,1,:));
        session_smoothed_tfs(:,:,end+1:end+size(img,3))=img;
    end
    session_mean_tfs=squeeze(mean(session_tfs,3));
    session_smoothed_mean_tfs=squeeze(mean(session_smoothed_tfs,3));
    
    mean_header=X(1);
    mean_header.fname=fullfile(epoch_dir, 'condition_Undefined.nii');
    spm_write_vol(mean_header, reshape(session_mean_tfs,size(session_mean_tfs,1),size(session_mean_tfs,2),1));
    
    mean_header=sX(1);
    mean_header.fname=fullfile(epoch_dir, 'scondition_Undefined.nii');
    spm_write_vol(mean_header, reshape(session_smoothed_mean_tfs,size(session_smoothed_mean_tfs,1),size(session_smoothed_mean_tfs,2),1));
end