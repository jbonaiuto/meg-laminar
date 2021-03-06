function extract_tfs(subj_info, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

subj_dir=fullfile('C:\meg_laminar\derivatives\spm12\', subj_info.subj_id);

dots_woi=[-3000 -500];
instr_woi=[-500 500];
resp_woi=[-1000 1000];

resp_channels={'MLC17','MLC25','MLC32','MLC42','MLC54','MLC63','MRC63',...
    'MLP57','MLP45','MLP35','MLP12','MLP23','MLC55','MZC04','MLP44',...
    'MLP34','MLP22','MLP11'};

dots_channels={'MLO53','MLO43','MLO32','MLO52','MLO31','MLO51','MLO41',...
    'MZO02','MZO03','MRO52','MRO42','MRO31','MRO53','MRO43','MRO32'};

instr_channels={'MLO53','MLO43','MLO32','MLO52','MLO31','MLO51','MLO41',...
    'MZO02','MZO03','MRO52','MRO42','MRO31','MRO53','MRO43','MRO32'};

sessions=[1:length(subj_info.sessions)];
if strcmp(subj_info.subj_id,'nc')
    sessions=[3];
end

for session_idx=1:length(sessions)
    session_num=sessions(session_idx);
    session_dir=fullfile(subj_dir, sprintf('ses-%02d',session_num));
    
    if exist(fullfile(session_dir, sprintf('rdots_tf_ffrcinstr_Tafdf%d.mat', session_num)),'file')==2
        spm_jobman('initcfg');
        matlabbatch={}; 
        clear jobs

        matlabbatch{1}.spm.meeg.images.convert2images.D = {fullfile(session_dir, sprintf('rdots_tf_ffrcinstr_Tafdf%d.mat', session_num))};
        matlabbatch{1}.spm.meeg.images.convert2images.mode = 'time x frequency';
        matlabbatch{1}.spm.meeg.images.convert2images.conditions = {};
        for idx=1:length(dots_channels)
            matlabbatch{1}.spm.meeg.images.convert2images.channels{idx}.chan = dots_channels{idx};
        end
        matlabbatch{1}.spm.meeg.images.convert2images.timewin = dots_woi;
        matlabbatch{1}.spm.meeg.images.convert2images.freqwin = [-Inf Inf];
        matlabbatch{1}.spm.meeg.images.convert2images.prefix = '';
        matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Convert2Images: M/EEG exported images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
        matlabbatch{2}.spm.spatial.smooth.dtype = 0;
        matlabbatch{2}.spm.spatial.smooth.im = 0;
        matlabbatch{2}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',matlabbatch);

        if exist(fullfile(session_dir,sprintf('rdots_tf_rcinstr_Tafdf%d',session_num)),'dir')==7
            rmdir(fullfile(session_dir,sprintf('rdots_tf_rcinstr_Tafdf%d',session_num)),'s');
        end
        delete(fullfile(session_dir, sprintf('rdots_tf_ffrcinstr_Tafdf%d.mat', session_num)));
        delete(fullfile(session_dir, sprintf('rdots_tf_ffrcinstr_Tafdf%d.dat', session_num)));
    end
end

for session_idx=1:length(sessions)
    session_num=sessions(session_idx);
    session_dir=fullfile(subj_dir, sprintf('ses-%02d',session_num));
    if exist(fullfile(session_dir, sprintf('rinstr_tf_ffrcinstr_Tafdf%d.mat', session_num)),'file')==2
        spm_jobman('initcfg');
        matlabbatch={}; 
        clear jobs

        matlabbatch{1}.spm.meeg.images.convert2images.D = {fullfile(session_dir, sprintf('rinstr_tf_ffrcinstr_Tafdf%d.mat', session_num))};
        matlabbatch{1}.spm.meeg.images.convert2images.mode = 'time x frequency';
        matlabbatch{1}.spm.meeg.images.convert2images.conditions = {};
        for idx=1:length(instr_channels)
            matlabbatch{1}.spm.meeg.images.convert2images.channels{idx}.chan = instr_channels{idx};
        end
        matlabbatch{1}.spm.meeg.images.convert2images.timewin = instr_woi;
        matlabbatch{1}.spm.meeg.images.convert2images.freqwin = [-Inf Inf];
        matlabbatch{1}.spm.meeg.images.convert2images.prefix = '';
        matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Convert2Images: M/EEG exported images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
        matlabbatch{2}.spm.spatial.smooth.dtype = 0;
        matlabbatch{2}.spm.spatial.smooth.im = 0;
        matlabbatch{2}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',matlabbatch);

        if exist(fullfile(session_dir,sprintf('rinstr_tf_rcinstr_Tafdf%d',session_num)),'dir')==7
            rmdir(fullfile(session_dir,sprintf('rinstr_tf_rcinstr_Tafdf%d',session_num)),'s');
        end
        delete(fullfile(session_dir, sprintf('rinstr_tf_ffrcinstr_Tafdf%d.mat', session_num)));
        delete(fullfile(session_dir, sprintf('rinstr_tf_ffrcinstr_Tafdf%d.dat', session_num)));
    end
end

for session_idx=1:length(sessions)
    session_num=sessions(session_idx);
    session_dir=fullfile(subj_dir, sprintf('ses-%02d',session_num));
    if exist(fullfile(session_dir, sprintf('rresp_tf_ffrcresp_Tafdf%d.mat', session_num)),'file')==2
        spm_jobman('initcfg');
        matlabbatch={}; 
        clear jobs

        matlabbatch{1}.spm.meeg.images.convert2images.D = {fullfile(session_dir, sprintf('rresp_tf_ffrcresp_Tafdf%d.mat', session_num))};
        matlabbatch{1}.spm.meeg.images.convert2images.mode = 'time x frequency';
        matlabbatch{1}.spm.meeg.images.convert2images.conditions = {};
        for idx=1:length(resp_channels)
            matlabbatch{1}.spm.meeg.images.convert2images.channels{idx}.chan = resp_channels{idx};
        end
        matlabbatch{1}.spm.meeg.images.convert2images.timewin = resp_woi;
        matlabbatch{1}.spm.meeg.images.convert2images.freqwin = [-Inf Inf];
        matlabbatch{1}.spm.meeg.images.convert2images.prefix = '';
        matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Convert2Images: M/EEG exported images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
        matlabbatch{2}.spm.spatial.smooth.dtype = 0;
        matlabbatch{2}.spm.spatial.smooth.im = 0;
        matlabbatch{2}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',matlabbatch);

        if exist(fullfile(session_dir,sprintf('rresp_tf_rcresp_Tafdf%d',session_num)),'dir')==7
            rmdir(fullfile(session_dir,sprintf('rresp_tf_rcresp_Tafdf%d',session_num)),'s');
        end
        delete(fullfile(session_dir, sprintf('rresp_tf_ffrcresp_Tafdf%d.mat', session_num)));
        delete(fullfile(session_dir, sprintf('rresp_tf_ffrcresp_Tafdf%d.dat', session_num)));
    end
end