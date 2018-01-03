function run_tfs(subj_info, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

subj_dir=fullfile('d:\pred_coding\derivatives\spm12', subj_info.subj_id);

dots_tf_woi=[-3500 -1000];
dots_baseline=[-3000 -2500];

instr_tf_woi=[-3500 1000];
instr_baseline=[-3000 -2500];

resp_tf_woi=[-1500 1500];
resp_baseline=[-3000 -2500];

for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-0%d',session_num));
    spm_jobman('initcfg');
    matlabbatch={}; 
    clear jobs

    matlabbatch{1}.spm.meeg.tf.tf.D = {fullfile(session_dir, sprintf('rcinstr_Tafdf%d.mat', session_num))};
    %matlabbatch{1}.spm.meeg.tf.tf.channels{1}.chan = dots_ch;
    matlabbatch{1}.spm.meeg.tf.tf.channels{1}.type = 'MEG';
    matlabbatch{1}.spm.meeg.tf.tf.frequencies = [2:100];
    matlabbatch{1}.spm.meeg.tf.tf.timewin = dots_tf_woi;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.ncycles = 7;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.timeres = 0;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.subsample = 1;
    matlabbatch{1}.spm.meeg.tf.tf.phase = 0;
    matlabbatch{1}.spm.meeg.tf.tf.prefix = 'tf_dots';
    matlabbatch{2}.spm.meeg.tf.rescale.D(1) = cfg_dep('Time-frequency analysis: M/EEG time-frequency power dataset', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dtfname'));
    matlabbatch{2}.spm.meeg.tf.rescale.method.LogR.baseline.timewin = dots_baseline;
    matlabbatch{2}.spm.meeg.tf.rescale.method.LogR.baseline.Db = [];
    matlabbatch{2}.spm.meeg.tf.rescale.prefix = 'r';
    spm_jobman('run',matlabbatch);
    
    delete(fullfile(session_dir, sprintf('tf_dotstf_rcinstr_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('tf_dotstf_rcinstr_Tafdf%d.dat', session_num)));
end

for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-0%d',session_num));
    spm_jobman('initcfg');
    matlabbatch={}; 
    clear jobs

    matlabbatch{1}.spm.meeg.tf.tf.D = {fullfile(session_dir, sprintf('rcinstr_Tafdf%d.mat', session_num))};
    %matlabbatch{1}.spm.meeg.tf.tf.channels{1}.chan = instr_ch;
    matlabbatch{1}.spm.meeg.tf.tf.channels{1}.type = 'MEG';
    matlabbatch{1}.spm.meeg.tf.tf.frequencies = [2:100];
    matlabbatch{1}.spm.meeg.tf.tf.timewin = instr_tf_woi;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.ncycles = 7;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.timeres = 0;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.subsample = 1;
    matlabbatch{1}.spm.meeg.tf.tf.phase = 0;
    matlabbatch{1}.spm.meeg.tf.tf.prefix = 'tf_instr';
    matlabbatch{2}.spm.meeg.tf.rescale.D(1) = cfg_dep('Time-frequency analysis: M/EEG time-frequency power dataset', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dtfname'));
    matlabbatch{2}.spm.meeg.tf.rescale.method.LogR.baseline.timewin = instr_baseline;
    matlabbatch{2}.spm.meeg.tf.rescale.method.LogR.baseline.Db = [];
    matlabbatch{2}.spm.meeg.tf.rescale.prefix = 'r';
    spm_jobman('run',matlabbatch);       
end

for session_num=1:length(subj_info.sessions)
    session_dir=fullfile(subj_dir, sprintf('ses-0%d',session_num));
    spm_jobman('initcfg');
    matlabbatch={}; 
    clear jobs

    matlabbatch{1}.spm.meeg.tf.tf.D = {fullfile(session_dir, sprintf('rcresp_Tafdf%d.mat', session_num))};
    %matlabbatch{1}.spm.meeg.tf.tf.channels{1}.chan = resp_ch;
    matlabbatch{1}.spm.meeg.tf.tf.channels{1}.type = 'MEG';
    matlabbatch{1}.spm.meeg.tf.tf.frequencies = [2:100];
    matlabbatch{1}.spm.meeg.tf.tf.timewin = resp_tf_woi;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.ncycles = 7;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.timeres = 0;
    matlabbatch{1}.spm.meeg.tf.tf.method.morlet.subsample = 1;
    matlabbatch{1}.spm.meeg.tf.tf.phase = 0;
    matlabbatch{1}.spm.meeg.tf.tf.prefix = '';
    matlabbatch{2}.spm.meeg.tf.rescale.D(1) = cfg_dep('Time-frequency analysis: M/EEG time-frequency power dataset', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dtfname'));
    matlabbatch{2}.spm.meeg.tf.rescale.method.LogR.baseline.timewin = resp_baseline;
    matlabbatch{2}.spm.meeg.tf.rescale.method.LogR.baseline.Db = {fullfile(session_dir, sprintf('tf_instrtf_rcinstr_Tafdf%d.mat', session_num))};
    matlabbatch{2}.spm.meeg.tf.rescale.prefix = 'r';
    spm_jobman('run',matlabbatch);
    
    delete(fullfile(session_dir, sprintf('tf_rcresp_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('tf_rcresp_Tafdf%d.dat', session_num)));
    delete(fullfile(session_dir, sprintf('tf_instrtf_rcinstr_Tafdf%d.mat', session_num)));
    delete(fullfile(session_dir, sprintf('tf_instrtf_rcinstr_Tafdf%d.dat', session_num)));
end
