function process_run_tfs(subj_info, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

subj_source_dir=fullfile('D:\meg_laminar\derivatives\spm12\', subj_info.subj_id);
subj_dest_dir=fullfile('C:\meg_laminar\derivatives\spm12\', subj_info.subj_id);

dots_tf_woi=[-3250 -250];
dots_baseline=[-3000 -2500];

instr_tf_woi=[-750 750];
instr_baseline=[-3000 -2500];

resp_tf_woi=[-1250 1250];
resp_baseline=[-3000 -2500];

sessions=[1:length(subj_info.sessions)];
if strcmp(subj_info.subj_id,'nc')
    sessions=[3];
end

for session_idx=1:length(sessions)
    session_num=sessions(session_idx);
    session_source_dir=fullfile(subj_source_dir, sprintf('ses-%02d',session_num));
    session_dest_dir=fullfile(subj_dest_dir, sprintf('ses-%02d',session_num));
    instr_fname=sprintf('rcinstr_Tafdf%d.mat', session_num);
    resp_fname=sprintf('rcresp_Tafdf%d.mat', session_num);
    
    spm_jobman('initcfg');
    
    % Copy file to foi_dir
    clear jobs
    matlabbatch={};
    matlabbatch{1}.spm.meeg.other.copy.D = {fullfile(session_source_dir, instr_fname)};
    matlabbatch{1}.spm.meeg.other.copy.outfile = fullfile(session_dest_dir, instr_fname);
    matlabbatch{2}.spm.meeg.other.copy.D = {fullfile(session_source_dir, resp_fname)};
    matlabbatch{2}.spm.meeg.other.copy.outfile = fullfile(session_dest_dir, resp_fname);
    spm_jobman('run',matlabbatch);
    
    % Relabel trials to all be same condition
    load(fullfile(session_dest_dir, instr_fname));
    D.condlist={''};
    for trial_idx=1:length(D.trials)
        D.trials(trial_idx).label='';
    end
    save(fullfile(session_dest_dir, instr_fname),'D');
    
    load(fullfile(session_dest_dir, resp_fname));
    D.condlist={''};
    for trial_idx=1:length(D.trials)
        D.trials(trial_idx).label='';
    end
    save(fullfile(session_dest_dir, resp_fname),'D');

    % Remove line noise
    clear jobs
    matlabbatch={};
    matlabbatch{1}.spm.meeg.preproc.filter.D = {fullfile(session_dest_dir, instr_fname)};
    matlabbatch{1}.spm.meeg.preproc.filter.type = 'butterworth';
    matlabbatch{1}.spm.meeg.preproc.filter.band = 'stop';
    matlabbatch{1}.spm.meeg.preproc.filter.freq = [49 51];
    matlabbatch{1}.spm.meeg.preproc.filter.dir = 'twopass';
    matlabbatch{1}.spm.meeg.preproc.filter.order = 5;
    matlabbatch{1}.spm.meeg.preproc.filter.prefix = 'f';
    matlabbatch{2}.spm.meeg.preproc.filter.D(1) = cfg_dep('Filter: Filtered Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
    matlabbatch{2}.spm.meeg.preproc.filter.type = 'butterworth';
    matlabbatch{2}.spm.meeg.preproc.filter.band = 'stop';
    matlabbatch{2}.spm.meeg.preproc.filter.freq = [99 101];
    matlabbatch{2}.spm.meeg.preproc.filter.dir = 'twopass';
    matlabbatch{2}.spm.meeg.preproc.filter.order = 5;
    matlabbatch{2}.spm.meeg.preproc.filter.prefix = 'f';
    spm_jobman('run',matlabbatch);
    clear jobs
    matlabbatch={};
    matlabbatch{1}.spm.meeg.preproc.filter.D = {fullfile(session_dest_dir, resp_fname)};
    matlabbatch{1}.spm.meeg.preproc.filter.type = 'butterworth';
    matlabbatch{1}.spm.meeg.preproc.filter.band = 'stop';
    matlabbatch{1}.spm.meeg.preproc.filter.freq = [49 51];
    matlabbatch{1}.spm.meeg.preproc.filter.dir = 'twopass';
    matlabbatch{1}.spm.meeg.preproc.filter.order = 5;
    matlabbatch{1}.spm.meeg.preproc.filter.prefix = 'f';
    matlabbatch{2}.spm.meeg.preproc.filter.D(1) = cfg_dep('Filter: Filtered Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
    matlabbatch{2}.spm.meeg.preproc.filter.type = 'butterworth';
    matlabbatch{2}.spm.meeg.preproc.filter.band = 'stop';
    matlabbatch{2}.spm.meeg.preproc.filter.freq = [99 101];
    matlabbatch{2}.spm.meeg.preproc.filter.dir = 'twopass';
    matlabbatch{2}.spm.meeg.preproc.filter.order = 5;
    matlabbatch{2}.spm.meeg.preproc.filter.prefix = 'f';
    spm_jobman('run',matlabbatch);
    
    % Dots
    session_dir=fullfile(subj_dest_dir, num2str(session_num));
    fname=sprintf('ffrcinstr_Tafdf%d.mat', session_num);
    
    spm_jobman('initcfg');    
    matlabbatch={}; 
    clear jobs

    matlabbatch{1}.spm.meeg.tf.tf.D = {fullfile(session_dir, fname)};
    matlabbatch{1}.spm.meeg.tf.tf.channels{1}.type = 'MEG';
    matlabbatch{1}.spm.meeg.tf.tf.frequencies = [2:100];
    matlabbatch{1}.spm.meeg.tf.tf.timewin = dots_tf_woi;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.ncycles = 7;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.timeres = 0;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.subsample = 1;
    matlabbatch{1}.spm.meeg.tf.tf.phase = 0;
    matlabbatch{1}.spm.meeg.tf.tf.prefix = 'dots_';
    matlabbatch{2}.spm.meeg.tf.rescale.D(1) = cfg_dep('Time-frequency analysis: M/EEG time-frequency power dataset', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dtfname'));
    matlabbatch{2}.spm.meeg.tf.rescale.method.Rel.baseline.timewin = dots_baseline;
    matlabbatch{2}.spm.meeg.tf.rescale.method.Rel.baseline.Db = [];
    matlabbatch{2}.spm.meeg.tf.rescale.prefix = 'r';
    spm_jobman('run',matlabbatch);
    
    
    % Instr
    session_dir=fullfile(subj_dest_dir, num2str(session_num));
    fname=sprintf('ffrcinstr_Tafdf%d.mat', session_num);
    
    spm_jobman('initcfg');
    matlabbatch={}; 
    clear jobs

    matlabbatch{1}.spm.meeg.tf.tf.D = {fullfile(session_dir, fname)};
    %matlabbatch{1}.spm.meeg.tf.tf.channels{1}.chan = instr_ch;
    matlabbatch{1}.spm.meeg.tf.tf.channels{1}.type = 'MEG';
    matlabbatch{1}.spm.meeg.tf.tf.frequencies = [2:100];
    matlabbatch{1}.spm.meeg.tf.tf.timewin = instr_tf_woi;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.ncycles = 7;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.timeres = 0;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.subsample = 1;
    matlabbatch{1}.spm.meeg.tf.tf.phase = 0;
    matlabbatch{1}.spm.meeg.tf.tf.prefix = 'instr_';
    matlabbatch{2}.spm.meeg.tf.rescale.D(1) = cfg_dep('Time-frequency analysis: M/EEG time-frequency power dataset', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dtfname'));
    matlabbatch{2}.spm.meeg.tf.rescale.method.Rel.baseline.timewin = instr_baseline;
    matlabbatch{2}.spm.meeg.tf.rescale.method.Rel.baseline.Db = {fullfile(session_dir, sprintf('dots_tf_ffrcinstr_Tafdf%d.mat', session_num))};
    matlabbatch{2}.spm.meeg.tf.rescale.prefix = 'r';
    spm_jobman('run',matlabbatch);       
    
    delete(fullfile(session_dir, sprintf('instr_tf_ffrcinstr_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('instr_tf_ffrcinstr_Tafdf%d.dat', session_num)));
    
    delete(fullfile(session_dir, sprintf('instr_tf_rcinstr_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('instr_tf_rcinstr_Tafdf%d.dat', session_num)));

    % Resp
    session_dir=fullfile(subj_dest_dir, num2str(session_num));
    fname=sprintf('ffrcresp_Tafdf%d.mat', session_num);
    
    spm_jobman('initcfg');
    matlabbatch={}; 
    clear jobs

    matlabbatch{1}.spm.meeg.tf.tf.D = {fullfile(session_dir, fname)};
    %matlabbatch{1}.spm.meeg.tf.tf.channels{1}.chan = resp_ch;
    matlabbatch{1}.spm.meeg.tf.tf.channels{1}.type = 'MEG';
    matlabbatch{1}.spm.meeg.tf.tf.frequencies = [2:100];
    matlabbatch{1}.spm.meeg.tf.tf.timewin = resp_tf_woi;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.ncycles = 7;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.timeres = 0;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.subsample = 1;
    matlabbatch{1}.spm.meeg.tf.tf.phase = 0;
    matlabbatch{1}.spm.meeg.tf.tf.prefix = 'resp_';
    matlabbatch{2}.spm.meeg.tf.rescale.D(1) = cfg_dep('Time-frequency analysis: M/EEG time-frequency power dataset', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dtfname'));
    matlabbatch{2}.spm.meeg.tf.rescale.method.Rel.baseline.timewin = resp_baseline;
    matlabbatch{2}.spm.meeg.tf.rescale.method.Rel.baseline.Db = {fullfile(session_dir, sprintf('dots_tf_ffrcinstr_Tafdf%d.mat', session_num))};
    matlabbatch{2}.spm.meeg.tf.rescale.prefix = 'r';
    spm_jobman('run',matlabbatch);
    
    delete(fullfile(session_dir, sprintf('rcresp_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('rcresp_Tafdf%d.dat', session_num)));
    delete(fullfile(session_dir, sprintf('frcresp_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('frcresp_Tafdf%d.dat', session_num)));
    delete(fullfile(session_dir, sprintf('ffrcresp_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('ffrcresp_Tafdf%d.dat', session_num)));
    
    delete(fullfile(session_dir, sprintf('rcinstr_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('rcinstr_Tafdf%d.dat', session_num)));
    delete(fullfile(session_dir, sprintf('frcinstr_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('frcinstr_Tafdf%d.dat', session_num)));
    delete(fullfile(session_dir, sprintf('ffrcinstr_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('ffrcinstr_Tafdf%d.dat', session_num)));
    
    delete(fullfile(session_dir, sprintf('resp_tf_rcresp_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('resp_tf_rcresp_Tafdf%d.dat', session_num)));
    delete(fullfile(session_dir, sprintf('dots_tf_rcinstr_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('dots_tf_rcinstr_Tafdf%d.dat', session_num)));
    
    delete(fullfile(session_dir, sprintf('resp_tf_ffrcresp_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('resp_tf_ffrcresp_Tafdf%d.dat', session_num)));
    delete(fullfile(session_dir, sprintf('dots_tf_ffrcinstr_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('dots_tf_ffrcinstr_Tafdf%d.dat', session_num)));
end
