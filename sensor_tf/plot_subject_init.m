function plot_subject_init(subj_info, zero_event, varargin)

defaults = struct('data_dir', '/data/pred_coding');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults', 'EEG');
spm_jobman('initcfg');

freq_bands=[7 13;15 35;60 90;2 100];
band_names={'alpha','beta','gamma','broadband'};
            
plot_analysis('scalp_freq',3);
plot_analysis('time_freq',4);

    function plot_analysis(type, units)
        analysis_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id, [type '_rtf_rc' zero_event '_Tafdf']);        
        dof=0;
        for i=1:length(subj_info.sessions)
            dof=dof+subj_info.sessions(i)*180;
        end
        V=spm_vol(fullfile(analysis_dir, 'spmT_0001.nii'));
        TMap = spm_read_vols(V);
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
                TMap,tthresh,fullfile(analysis_dir,'t_test_positive.png'), freq_bands, event_times);            
            
        elseif units==3
            coords=V.mat*[[1:max_dim]' [1:max_dim]' [1:max_dim]' ones(max_dim,1)]';
            
            for i=1:size(freq_bands,1)
                plot_scalp_frequency(coords(2,1:V.dim(2)),...
                    coords(1,1:V.dim(1)),coords(3,1:V.dim(3)),TMap,...
                    tthresh,freq_bands(i,:),...
                    fullfile(analysis_dir,sprintf('t_test_positive_%s.png',band_names{i})))                               
            end
        end
                
        V=spm_vol(fullfile(analysis_dir, 'spmT_0002.nii'));
        TMap = spm_read_vols(V);
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
                TMap,tthresh,fullfile(analysis_dir,'t_test_negative.png'), freq_bands, event_times);            
            
        elseif units==3
            coords=V.mat*[[1:max_dim]' [1:max_dim]' [1:max_dim]' ones(max_dim,1)]';
            
            for i=1:size(freq_bands,1)
                plot_scalp_frequency(coords(2,1:V.dim(2)),...
                    coords(1,1:V.dim(1)),coords(3,1:V.dim(3)),TMap,...
                    tthresh,freq_bands(i,:),...
                    fullfile(analysis_dir,sprintf('t_test_negative_%s.png',band_names{i})))                
            end
        end
        
    end

    function plot_time_frequency(times,freqs, x, threshold, filename, freq_bands, event_times)
        f=figure();
        imagesc(times,freqs,x);
        set(gca,'clim',[threshold max(x(:))]);
        set(gca,'ydir','normal');
        hold on;
        for i=1:length(event_times)
            plot([event_times(i) event_times(i)],[freqs(1) freqs(end)],'w--');
        end
        for i=1:size(freq_bands,1)-1
            plot([times(1) times(end)],[freq_bands(i,1) freq_bands(i,1)],'w--');
            plot([times(1) times(end)],[freq_bands(i,2) freq_bands(i,2)],'w--');
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
        imagesc(y,freqs(freq_idx),yf');
        set(gca,'clim',[threshold max([threshold+1 max(yf(:))])]);
        set(gca,'ydir','normal');
        title(sprintf('%d-%dHz',freq_band(1),freq_band(2)));
        xlabel('Posterior - Anterior');
        ylabel('Frequency (Hz)');
        colorbar();
        
        subplot(2,2,2);
        imagesc(y,freqs(freq_idx),xf');
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
