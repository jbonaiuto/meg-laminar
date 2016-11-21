function all_subjects_sensor_tf_init(subjects, zero_event, varargin)

defaults = struct('data_dir', '/data/pred_coding');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults', 'EEG');
spm_jobman('initcfg');

for subj_idx=1:length(subjects)
    plot_subject_init(subjects(subj_idx),zero_event,'data_dir',params.data_dir);
end

conditions={'congruent-low','congruent-med','congruent-high',...
    'incongruent-low','incongruent-med','incongruent-high'};

run_analysis('scalp_freq',3);
%run_analysis('scalp_time',2);
run_analysis('time_freq',4);

    function run_analysis(type, units)
        curr_dir=pwd;
        spm_jobman('initcfg');
        clear jobs;
        matlabbatch={};
        batch_idx=1;
        
        analysis_dir=fullfile(params.data_dir, 'analysis', [type '_rtf_rc' zero_event '_Tafdf']);
        delete(fullfile(analysis_dir,'*'));
        
        matlabbatch{batch_idx}.spm.stats.factorial_design.dir = {analysis_dir};
        matlabbatch{batch_idx}.spm.stats.factorial_design.des.fd.fact(1).name = 'direction';
        matlabbatch{batch_idx}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
        matlabbatch{batch_idx}.spm.stats.factorial_design.des.fd.fact(1).dept = 0;
        matlabbatch{batch_idx}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
        matlabbatch{batch_idx}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
        matlabbatch{batch_idx}.spm.stats.factorial_design.des.fd.icell(1).levels = 1;
        %%
        matlabbatch{batch_idx}.spm.stats.factorial_design.des.fd.icell(1).scans = {};
        for subj_idx=1:length(subjects)
            subj_info=subjects(subj_idx);
            matlabbatch{batch_idx}.spm.stats.factorial_design.des.fd.icell(1).scans{end+1,1}=fullfile(params.data_dir, 'analysis', ...
                subj_info.subj_id, [type '_rtf_rc' zero_event '_Tafdf'],'con_0001.nii,1');
        end
        %%
        matlabbatch{batch_idx}.spm.stats.factorial_design.des.fd.icell(2).levels = 2;
        %%
        matlabbatch{batch_idx}.spm.stats.factorial_design.des.fd.icell(2).scans = {};
        for subj_idx=1:length(subjects)
            subj_info=subjects(subj_idx);
            matlabbatch{batch_idx}.spm.stats.factorial_design.des.fd.icell(2).scans{end+1,1}=fullfile(params.data_dir, 'analysis', ...
                subj_info.subj_id, [type '_rtf_rc' zero_event '_Tafdf'],'con_0002.nii,1');
        end
        %%
        matlabbatch{batch_idx}.spm.stats.factorial_design.des.fd.contrasts = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{batch_idx}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{batch_idx}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{batch_idx}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.glonorm = 1;
        batch_idx=batch_idx+1;
        
        matlabbatch{batch_idx}.spm.stats.fmri_est.spmmat(1) = {fullfile(analysis_dir,'SPM.mat')};
        matlabbatch{batch_idx}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{batch_idx}.spm.stats.fmri_est.method.Classical = 1;
        batch_idx=batch_idx+1;
        
        matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = {fullfile(analysis_dir,'SPM.mat')};
        matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.name = 'average';
        matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.weights = [1 0;0 1];
        matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
        matlabbatch{batch_idx}.spm.stats.con.delete = 1;
        batch_idx=batch_idx+1;
        
        % matlabbatch{batch_idx}.spm.stats.results.spmmat(1) = {fullfile(analysis_dir,'SPM.mat')};
        % matlabbatch{batch_idx}.spm.stats.results.conspec.titlestr = '';
        % matlabbatch{batch_idx}.spm.stats.results.conspec.contrasts = Inf;
        % matlabbatch{batch_idx}.spm.stats.results.conspec.threshdesc = 'none';
        % matlabbatch{batch_idx}.spm.stats.results.conspec.thresh = 0.05;
        % matlabbatch{batch_idx}.spm.stats.results.conspec.extent = 0;
        % matlabbatch{batch_idx}.spm.stats.results.conspec.mask.none = 1;
        % matlabbatch{batch_idx}.spm.stats.results.units = units;
        % matlabbatch{batch_idx}.spm.stats.results.print = 'png';
        % matlabbatch{batch_idx}.spm.stats.results.write.none = 1;
        % batch_idx=batch_idx+1;
        
        spm_jobman('run', matlabbatch);
        
        V=spm_vol(fullfile(analysis_dir, 'spmF_0001.nii'));
        FMap = spm_read_vols(V);
        dof1=2;
        dof2=length(subjects)*2-2;
        fthresh = finv(.95,dof1,dof2);  % This is the corresponding t-threshold I'll need to extract all voxels whose value is greater than this.
        max_dim=max(V.dim);
        
        if units==4
            coords=V.mat*[[1:max_dim]' [1:max_dim]' ones(max_dim,1) ones(max_dim,1)]';
            if strcmp(zero_event,'dots')
                coords(2,:)=coords(2,:)+2500;
            end
            event_times=[0];
            if strcmp(zero_event,'instr')
                event_times(end+1)=-2500;
            end
            plot_time_frequency(coords(2,1:V.dim(2)),coords(1,1:V.dim(1)),...
                FMap,fthresh,fullfile(analysis_dir,'f_test.png'), event_times);            
        elseif units==3
            coords=V.mat*[[1:max_dim]' [1:max_dim]' [1:max_dim]' ones(max_dim,1)]';
            
            freq_bands=[15 35;60 90;2 100];
            band_names={'beta','gamma','broadband'};
            
            for i=1:size(freq_bands,1)
                plot_scalp_frequency(coords(2,1:V.dim(2)),...
                    coords(1,1:V.dim(1)),coords(3,1:V.dim(3)),FMap,...
                    fthresh,freq_bands(i,:),...
                    fullfile(analysis_dir,sprintf('f_test_%s.png',band_names{i})))                
            end
        end
        
        clear jobs;
        matlabbatch={};
        batch_idx=1;
        
        analysis_dir=fullfile(params.data_dir, 'analysis', [type '_rtf_rc' zero_event '_Tafdf_positive']);
        delete(fullfile(analysis_dir,'*'));
        
        matlabbatch{batch_idx}.spm.stats.factorial_design.dir = {analysis_dir};
        matlabbatch{batch_idx}.spm.stats.factorial_design.des.t1.scans = {};
        for subj_idx=1:length(subjects)
            subj_info=subjects(subj_idx);
            matlabbatch{batch_idx}.spm.stats.factorial_design.des.t1.scans{end+1,1}=fullfile(params.data_dir, 'analysis', ...
                subj_info.subj_id, [type '_rtf_rc' zero_event '_Tafdf'],'con_0001.nii,1');
        end
        matlabbatch{batch_idx}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{batch_idx}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{batch_idx}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{batch_idx}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.glonorm = 1;
        batch_idx=batch_idx+1;
        
        matlabbatch{batch_idx}.spm.stats.fmri_est.spmmat(1) = {fullfile(analysis_dir,'SPM.mat')};
        matlabbatch{batch_idx}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{batch_idx}.spm.stats.fmri_est.method.Classical = 1;
        batch_idx=batch_idx+1;
        
        matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = {fullfile(analysis_dir,'SPM.mat')};
        matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.name = 'positive';
        matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{batch_idx}.spm.stats.con.delete = 1;
        batch_idx=batch_idx+1;
        
        % matlabbatch{batch_idx}.spm.stats.results.spmmat(1) = {fullfile(analysis_dir,'SPM.mat')};
        % matlabbatch{batch_idx}.spm.stats.results.conspec.titlestr = '';
        % matlabbatch{batch_idx}.spm.stats.results.conspec.contrasts = Inf;
        % matlabbatch{batch_idx}.spm.stats.results.conspec.threshdesc = 'none';
        % matlabbatch{batch_idx}.spm.stats.results.conspec.thresh = 0.05;
        % matlabbatch{batch_idx}.spm.stats.results.conspec.extent = 0;
        % matlabbatch{batch_idx}.spm.stats.results.conspec.mask.none = 1;
        % matlabbatch{batch_idx}.spm.stats.results.units = units;
        % matlabbatch{batch_idx}.spm.stats.results.print = 'png';
        % matlabbatch{batch_idx}.spm.stats.results.write.none = 1;
        % batch_idx=batch_idx+1;
        
        spm_jobman('run', matlabbatch);
        
        V=spm_vol(fullfile(analysis_dir, 'spmT_0001.nii'));
        TMap = spm_read_vols(V);
        dof=length(subjects)-1;
        tthresh = tinv(.95,dof);  % This is the corresponding t-threshold I'll need to extract all voxels whose value is greater than this.
        max_dim=max(V.dim);
        
        if units==4
            coords=V.mat*[[1:max_dim]' [1:max_dim]' ones(max_dim,1) ones(max_dim,1)]';
            if strcmp(zero_event,'dots')
                coords(2,:)=coords(2,:)+2500;
            end
            event_times=[0];
            if strcmp(zero_event,'instr')
                event_times(end+1)=-2500;
            end
            plot_time_frequency(coords(2,1:V.dim(2)),coords(1,1:V.dim(1)),...
                TMap,tthresh,fullfile(analysis_dir,'t_test_positive.png'), event_times);            
            
        elseif units==3
            coords=V.mat*[[1:max_dim]' [1:max_dim]' [1:max_dim]' ones(max_dim,1)]';
            
            freq_bands=[15 35;60 90;2 100];
            band_names={'beta','gamma','broadband'};
            
            for i=1:size(freq_bands,1)
                plot_scalp_frequency(coords(2,1:V.dim(2)),...
                    coords(1,1:V.dim(1)),coords(3,1:V.dim(3)),TMap,...
                    tthresh,freq_bands(i,:),...
                    fullfile(analysis_dir,sprintf('t_test_positive_%s.png',band_names{i})))                               
            end
        end
        
        clear jobs;
        matlabbatch={};
        batch_idx=1;
        
        analysis_dir=fullfile(params.data_dir, 'analysis', [type '_rtf_rc' zero_event '_Tafdf_negative']);
        delete(fullfile(analysis_dir,'*'));
        
        matlabbatch{batch_idx}.spm.stats.factorial_design.dir = {analysis_dir};
        matlabbatch{batch_idx}.spm.stats.factorial_design.des.t1.scans = {};
        for subj_idx=1:length(subjects)
            subj_info=subjects(subj_idx);
            matlabbatch{batch_idx}.spm.stats.factorial_design.des.t1.scans{end+1,1}=fullfile(params.data_dir, 'analysis', ...
                subj_info.subj_id, [type '_rtf_rc' zero_event '_Tafdf'],'con_0002.nii,1');
        end
        matlabbatch{batch_idx}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{batch_idx}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{batch_idx}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{batch_idx}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.glonorm = 1;
        batch_idx=batch_idx+1;
        
        matlabbatch{batch_idx}.spm.stats.fmri_est.spmmat(1) = {fullfile(analysis_dir,'SPM.mat')};
        matlabbatch{batch_idx}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{batch_idx}.spm.stats.fmri_est.method.Classical = 1;
        batch_idx=batch_idx+1;
        
        matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = {fullfile(analysis_dir,'SPM.mat')};
        matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.name = 'negative';
        matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{batch_idx}.spm.stats.con.delete = 1;
        batch_idx=batch_idx+1;
        
        % matlabbatch{batch_idx}.spm.stats.results.spmmat(1) = {fullfile(analysis_dir,'SPM.mat')};
        % matlabbatch{batch_idx}.spm.stats.results.conspec.titlestr = '';
        % matlabbatch{batch_idx}.spm.stats.results.conspec.contrasts = Inf;
        % matlabbatch{batch_idx}.spm.stats.results.conspec.threshdesc = 'none';
        % matlabbatch{batch_idx}.spm.stats.results.conspec.thresh = 0.05;
        % matlabbatch{batch_idx}.spm.stats.results.conspec.extent = 0;
        % matlabbatch{batch_idx}.spm.stats.results.conspec.mask.none = 1;
        % matlabbatch{batch_idx}.spm.stats.results.units = units;
        % matlabbatch{batch_idx}.spm.stats.results.print = 'png';
        % matlabbatch{batch_idx}.spm.stats.results.write.none = 1;
        % batch_idx=batch_idx+1;
        
        spm_jobman('run', matlabbatch);
        
        V=spm_vol(fullfile(analysis_dir, 'spmT_0001.nii'));
        TMap = spm_read_vols(V);
        dof=length(subjects)-1;
        tthresh = tinv(.95,dof);  % This is the corresponding t-threshold I'll need to extract all voxels whose value is greater than this.
        max_dim=max(V.dim);
        if units==4
            coords=V.mat*[[1:max_dim]' [1:max_dim]' ones(max_dim,1) ones(max_dim,1)]';
            if strcmp(zero_event,'dots')
                coords(2,:)=coords(2,:)+2500;
            end
            event_times=[0];
            if strcmp(zero_event,'instr')
                event_times(end+1)=-2500;
            end
            plot_time_frequency(coords(2,1:V.dim(2)),coords(1,1:V.dim(1)),...
                TMap,tthresh,fullfile(analysis_dir,'t_test_negative.png'), event_times);            
            
        elseif units==3
            coords=V.mat*[[1:max_dim]' [1:max_dim]' [1:max_dim]' ones(max_dim,1)]';
            freq_bands=[15 35;60 90;2 100];
            band_names={'beta','gamma','broadband'};
            
            for i=1:size(freq_bands,1)
                plot_scalp_frequency(coords(2,1:V.dim(2)),...
                    coords(1,1:V.dim(1)),coords(3,1:V.dim(3)),TMap,...
                    tthresh,freq_bands(i,:),...
                    fullfile(analysis_dir,sprintf('t_test_negative_%s.png',band_names{i})))                
            end
        end
        
        cd(curr_dir);
    end

    function plot_time_frequency(times,freqs, x, threshold, filename, event_times)
        f=figure();
        imagesc(times,freqs,x);
        set(gca,'clim',[threshold max(x(:))]);
        set(gca,'ydir','normal');
        hold on;
        for i=1:length(event_times)
            plot([event_times(i) event_times(i)],[freqs(1) freqs(end)],'w--');
        end
        ylabel('Frequency (Hz)');
        xlabel('Time (ms)');
        colorbar();
        saveas(f,filename);
    end

    function plot_scalp_frequency(x,y,freqs,sf,threshold,freq_band,filename)
        freq_idx=intersect(find(freqs>=freq_band(1)),find(freqs<=freq_band(2)));
        yf=squeeze(max(sf(:,:,freq_idx),[],1));
        xf=squeeze(max(sf(:,:,freq_idx),[],2));
        xy=squeeze(max(sf(:,:,freq_idx),[],3));
        
        f=figure('position', [0, 0, 1200, 600]);
        subplot(2,2,1);
        imagesc(y,freqs,yf');
        set(gca,'clim',[threshold max([threshold+1 max(yf(:))])]);
        set(gca,'ydir','normal');
        title(sprintf('%d-%dHz',freq_band(1),freq_band(2)));
        xlabel('Posterior - Anterior');
        ylabel('Frequency (Hz)');
        colorbar();
        
        subplot(2,2,2);
        imagesc(y,freqs,xf');
        set(gca,'clim',[threshold max([threshold+1 max(xf(:))])]);
        set(gca,'ydir','normal');
        title(sprintf('%d-%dHz',freq_band(1),freq_band(2)));
        xlabel('Right - Left');
        ylabel('Frequency (Hz)');
        colorbar();
        
        subplot(2,2,3);
        imagesc(x,y,xy);
        set(gca,'clim',[threshold max([threshold+1 max(xy(:))])]);
        %set(gca,'ydir','normal');
        title(sprintf('%d-%dHz',freq_band(1),freq_band(2)));
        xlabel('Posterior - Anterior');
        ylabel('Right - Left');
        colorbar();
        saveas(f,filename);
    end
end
