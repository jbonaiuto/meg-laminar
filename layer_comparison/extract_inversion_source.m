function extract_inversion_source(subj_info, session_num, zero_event, foi, wois, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'inv_type', 'EBB', 'patch_size',0.4, 'parallel', true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

grey_coreg_dir=fullfile(params.data_dir,'analysis', subj_info.subj_id, num2str(session_num), 'grey_coreg');
foi_dir=fullfile(grey_coreg_dir, params.inv_type, ['p' num2str(params.patch_size)], zero_event, ['f' num2str(foi(1)) '_' num2str(foi(2))]);

coreg_file_name=fullfile(foi_dir, sprintf('r%s_%d.mat', subj_info.subj_id, session_num));

spm('defaults', 'EEG');

NEW=0;
if NEW
    D=spm_eeg_load(coreg_file_name);
    goodchans=D.indchantype('MEGGRAD','good');
    Dgood=squeeze(D(goodchans,:,:));
    M=D.inv{1}.inverse.M;
    U=D.inv{1}.inverse.U{1};
    MU=M*U;

    n_combined_vertices=size(M,1);
    n_vertices=round(n_combined_vertices/2);

    times=D.inv{1}.inverse.pst;
    ntrials=size(Dgood,3);
    for w=1:size(wois,1)
        woi_dir=fullfile(foi_dir, ['t' num2str(wois(w,1)) '_' num2str(wois(w,2))]);
        if exist(woi_dir,'dir')~=7
            mkdir(woi_dir);
        end
    end       
    woi_trial_vals=zeros(n_combined_vertices,ntrials);
    if params.parallel
        parfor t=1:ntrials
            d1=squeeze(Dgood(:,:,t));
            Dtrial=MU*d1;
            for w=1:size(wois,1)
                woi_idx=intersect(find(times>=wois(w,1)),find(times<=wois(w,2)));
                woi_trial_vals(:,t)=sum(Dtrial(:,woi_idx).^2,2);
            end
        end
    else
        for t=1:ntrials
            d1=squeeze(Dgood(:,:,t));
            Dtrial=MU*d1;
            for w=1:size(wois,1)
                woi_idx=intersect(find(times>=wois(w,1)),find(times<=wois(w,2)));
                woi_trial_vals(:,t)=sum(Dtrial(:,woi_idx).^2,2);
            end
        end
    end
    for w=1:size(wois,1)
        woi_dir=fullfile(foi_dir, ['t' num2str(wois(w,1)) '_' num2str(wois(w,2))]);
        for t=1:ntrials
            file_prefix=sprintf('pial_r%s_%d_1_t%d_%d_f%d_%d_1_%d',subj_info.subj_id,session_num,woi(1),woi(2),foi(1),foi(2),t);
            c=file_array(fullfile(woi_dir,sprintf('%s.dat',file_prefix)),[n_vertices 1],'FLOAT32-LE',0,1,0);
            c(:)=woi_trial_vals(n_vertices+1:end,t);
            pial_t=gifti;
            pial_t.cdata=c;
            save(pial_t, fullfile(woi_dir,sprintf('%s.gii',file_prefix)));

            file_prefix=sprintf('white_r%s_%d_1_t%d_%d_f%d_%d_1_%d',subj_info.subj_id,session_num,woi(1),woi(2),foi(1),foi(2),t);
            c=file_array(fullfile(woi_dir,sprintf('%s.dat',file_prefix)),[n_vertices 1],'FLOAT32-LE',0,1,0);
            c(:)=woi_trial_vals(1:n_vertices,t);
            white_t=gifti;
            white_t.cdata=c;
            save(white_t, fullfile(woi_dir,sprintf('%s.gii',file_prefix)));
        end

    end
else
    for w=1:size(wois,1)
        woi=wois(w,:);
        % Extract source at woi1 and woi2
        spm_jobman('initcfg'); 
        clear jobs
        matlabbatch={};
        batch_idx=1;

        matlabbatch{batch_idx}.spm.meeg.source.results.D = {coreg_file_name};
        matlabbatch{batch_idx}.spm.meeg.source.results.val = 1;
        matlabbatch{batch_idx}.spm.meeg.source.results.woi = woi;
        matlabbatch{batch_idx}.spm.meeg.source.results.foi = foi;
        matlabbatch{batch_idx}.spm.meeg.source.results.ctype = 'trials';
        matlabbatch{batch_idx}.spm.meeg.source.results.space = 0;
        matlabbatch{batch_idx}.spm.meeg.source.results.format = 'mesh';
        matlabbatch{batch_idx}.spm.meeg.source.results.smoothing = 8;
        batch_idx=batch_idx+1;

        spm_jobman('run',matlabbatch);

        % Move files to subdirectories
        woi_dir=fullfile(foi_dir, ['t' num2str(woi(1)) '_' num2str(woi(2))]);
        if exist(woi_dir,'dir')~=7
            mkdir(woi_dir);
        end
        movefile(fullfile(foi_dir, ['r' subj_info.subj_id '_' num2str(session_num) '_1_t' num2str(woi(1)) '_' num2str(woi(2)) '_f' num2str(foi(1)) '_' num2str(foi(2)) '_*']), woi_dir);

        % Split pial and grey sources
        split_inversion_results(subj_info, grey_coreg_dir, zero_event, foi, woi, 'patch_size', params.patch_size, 'inv_type', params.inv_type);
    end

end



