function results=plot_subject_tfs(subj_info, epoch_name, zero_evt, varargin)

defaults = struct('plot',true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

subj_dir=fullfile('C:\meg_laminar\derivatives\spm12\', subj_info.subj_id);

session_tfs=[];
sessions=[1:length(subj_info.sessions)];
if strcmp(subj_info.subj_id,'nc')
    sessions=[3];
end
for session_idx=1:length(sessions)
    session_num=sessions(session_idx);
    session_dir=fullfile(subj_dir, sprintf('ses-%02d',session_num));
    tf_dir=fullfile(session_dir,sprintf('r%s_tf_ffrc%s_Tafdf%d',epoch_name,zero_evt,session_num));
    if exist(tf_dir,'dir')~=7
        tf_dir=fullfile(session_dir,sprintf('r%s_tf_rc%s_Tafdf%d',epoch_name,zero_evt,session_num));
    end
    
    X=spm_vol(fullfile(tf_dir, 'condition_Undefined.nii'));
    max_dim=max(X(1).dim);
    coords=X(1).mat*[[1:max_dim]' [1:max_dim]' ones(max_dim,1) ones(max_dim,1)]';
    
    freq_idx=intersect(find(coords(1,1:X(1).dim(1))>=10),find(coords(1,1:X(1).dim(1))<=40));
    %freq_idx=[1:X(1).dim(1)];
    results.times=coords(2,1:X(1).dim(2));
    results.freqs=coords(1,freq_idx);
        
    img=spm_read_vols(X);
    img=squeeze(img(freq_idx,:,1,:));
    session_tfs(:,:,end+1:end+size(img,3))=10.*log10(img./100+1);
end
results.session_mean_tfs=squeeze(mean(session_tfs,3));

if params.plot
    fig=figure();
    imagesc(results.times, results.freqs, results.session_mean_tfs);
    set(gca,'ydir','normal');
    colorbar();
    xlabel('Time');
    ylabel('Freq');
end