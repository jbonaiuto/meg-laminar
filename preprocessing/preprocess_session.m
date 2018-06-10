function preprocess_session(subj_info, session_num, varargin)

defaults = struct('data_dir', 'd:\meg_laminar\', ...
    'convert',true, ...
    'filter', true, 'remove_blinks', true, 'epoch', true, 'delete_last_resp', false, ...
    'highpass_freq',2.0,'downsample',250,'lowpass_freq',100.0,'instr_epoch',[-3500 1500],'resp_epoch',[-2000 2000], ...
    'blink_channel','MLT31','manual_reject',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Add path to report generation code
addpath C:/users/jbonai/Documents/MATLAB/m2html/
    
% Directory containing session data
session_dir=fullfile(params.data_dir,subj_info.subj_id, sprintf('ses-0%d',session_num));
% Directory to put results
analysis_dir=fullfile(params.data_dir,'derivatives/spm12',subj_info.subj_id,sprintf('ses-0%d',session_num));
if exist(analysis_dir,'dir')~=7
    mkdir(analysis_dir);
end

spm('defaults','eeg');
spm_jobman('initcfg');

instr_run_files={};
resp_run_files={};

for run_num=1:subj_info.sessions(session_num)
    
    % Directory containing run data
    run_code=[subj_info.subj_id '_JamesBonaiuto_' subj_info.scan_date{session_num} '_0' num2str(run_num)];
    run_dir=fullfile(session_dir, 'meg', [run_code '.ds']);

    % File containing stimulus information for each trial
    stim_file=fullfile(session_dir, 'behavior', ['stim_' num2str(run_num) '.mat']);
    % File containing behavioral data from each trial
    data_file=fullfile(session_dir, 'behavior', ['data_' num2str(run_num) '.mat']);

    spm_filename=sprintf('%d-%d.mat',session_num,run_num);
  
    if params.convert
        clear jobs
        batch_idx=1;
        matlabbatch={};

        % Convert to SPM format
        matlabbatch{batch_idx}.spm.meeg.convert.dataset = {fullfile(run_dir, [run_code '.meg4'])};
        matlabbatch{batch_idx}.spm.meeg.convert.mode.continuous.readall = 1;
        matlabbatch{batch_idx}.spm.meeg.convert.channels{batch_idx}.all = 'all';
        matlabbatch{batch_idx}.spm.meeg.convert.outfile = fullfile(analysis_dir, spm_filename);
        matlabbatch{batch_idx}.spm.meeg.convert.eventpadding = 0;
        matlabbatch{batch_idx}.spm.meeg.convert.blocksize = 3276800;
        matlabbatch{batch_idx}.spm.meeg.convert.checkboundary = 1;
        matlabbatch{batch_idx}.spm.meeg.convert.saveorigheader = 0;
        matlabbatch{batch_idx}.spm.meeg.convert.inputformat = 'autodetect';
        spm_jobman('run',matlabbatch);

        plot_fiducial_coils(fullfile(analysis_dir, spm_filename), 'lim_last_event',false,'lim_jump',false);
        
        plot_fiducial_coils(fullfile(analysis_dir, spm_filename), 'lim_last_event',true,'lim_jump',false);
        
        plot_fiducial_coils(fullfile(analysis_dir, spm_filename), 'lim_last_event', false, 'lim_jump', true);
        
        % Correct event timings, adjust initial events
        [num_diode_onsets num_dots_evts num_instr_evts num_resp_evts]=adjust_run_events(subj_info, session_num, run_num, 'data_dir', params.data_dir, 'plot_diode',true,...
            'delete_last_resp',params.delete_last_resp,'delete_no_resp',true);        

    end
    
    if params.filter
        clear jobs
        matlabbatch={};
        batch_idx=1;

        % Highpass filter
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.D = {fullfile(analysis_dir, spm_filename)};
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.type = 'butterworth';
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.band = 'high';
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.freq = params.highpass_freq;
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.dir = 'twopass';
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.order = 5;
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.prefix = fullfile(analysis_dir,'f');
        batch_idx=batch_idx+1;
        
        % Downsample
        matlabbatch{batch_idx}.spm.meeg.preproc.downsample.D = {fullfile(analysis_dir, sprintf('f%s',spm_filename))};
        matlabbatch{batch_idx}.spm.meeg.preproc.downsample.fsample_new = params.downsample;
        matlabbatch{batch_idx}.spm.meeg.preproc.downsample.prefix = fullfile(analysis_dir,'d');
        batch_idx=batch_idx+1;
        
        % Lowpass filter
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.D = {fullfile(analysis_dir, sprintf('df%s',spm_filename))};
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.type = 'butterworth';
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.band = 'low';
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.freq = params.lowpass_freq;
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.dir = 'twopass';
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.order = 5;
        matlabbatch{batch_idx}.spm.meeg.preproc.filter.prefix = fullfile(analysis_dir,'f');
        batch_idx=batch_idx+1;
        spm_jobman('run',matlabbatch);
    end
    
    if params.remove_blinks
        % Remove blinks
        clear jobs
        matlabbatch={};
        batch_idx=1;

        % Coregister with MRI
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.D = {fullfile(analysis_dir, sprintf('fdf%s',spm_filename))};
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.val = 1;
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.comment = '';
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {fullfile(params.data_dir,'mri',[subj_info.subj_id subj_info.birth_date], [subj_info.headcast_t1 ',1'])};
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshres = 2;
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = subj_info.nas;
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = subj_info.lpa;
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = subj_info.rpa;
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
        matlabbatch{batch_idx}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
        batch_idx=batch_idx+1;

        matlabbatch{batch_idx}.spm.meeg.preproc.artefact.D = {fullfile(analysis_dir, sprintf('fdf%s',spm_filename))};
        matlabbatch{batch_idx}.spm.meeg.preproc.artefact.mode = 'mark';
        matlabbatch{batch_idx}.spm.meeg.preproc.artefact.badchanthresh = 0.2;
        matlabbatch{batch_idx}.spm.meeg.preproc.artefact.append = true;
        matlabbatch{batch_idx}.spm.meeg.preproc.artefact.methods.channels{1}.chan = params.blink_channel;
        matlabbatch{batch_idx}.spm.meeg.preproc.artefact.methods.fun.eyeblink.threshold = 4;
        matlabbatch{batch_idx}.spm.meeg.preproc.artefact.methods.fun.eyeblink.excwin = 0;
        matlabbatch{batch_idx}.spm.meeg.preproc.artefact.prefix = 'a';
        batch_idx=batch_idx+1;

        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.D = {fullfile(analysis_dir, sprintf('afdf%s',spm_filename))};
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.timewin = [-500 500];
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.conditionlabel = 'Eyeblink';
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.eventtype = 'artefact_eyeblink';
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.eventvalue = params.blink_channel;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.trlshift = 0;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.bc = 0;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.eventpadding = 0;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.prefix = 'eyeblink';
        batch_idx=batch_idx+1;

        if ~length(subj_info.blink_ncomps)
            matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.D = {fullfile(analysis_dir, sprintf('eyeblinkafdf%s',spm_filename))};
            matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.mode{1}.svd.timewin = [-Inf Inf];
            matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.mode{1}.svd.threshold = NaN;
            matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.mode{1}.svd.ncomp = 4;
            batch_idx=batch_idx+1;

            spm_jobman('run',matlabbatch);

            clear jobs
            matlabbatch={};
            batch_idx=1;

            ncomp=input('Number of components: ');
            matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.D = {fullfile(analysis_dir, sprintf('eyeblinkafdf%s',spm_filename))};
            matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.mode{1}.clear=1;
            batch_idx=batch_idx+1;
        else
            ncomp=subj_info.blink_ncomps(session_num,run_num);
        end

        matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.D = {fullfile(analysis_dir, sprintf('eyeblinkafdf%s',spm_filename))};
        matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.mode{1}.svd.timewin = [-Inf Inf];
        matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.mode{1}.svd.threshold = NaN;
        matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.mode{1}.svd.ncomp = ncomp;
        batch_idx=batch_idx+1;

        matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.D = {fullfile(analysis_dir, sprintf('afdf%s',spm_filename))};
        matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.mode{1}.spmeeg.conffile = {fullfile(analysis_dir, sprintf('eyeblinkafdf%s',spm_filename))};
        batch_idx=batch_idx+1;

        matlabbatch{batch_idx}.spm.meeg.preproc.correct.D = {fullfile(analysis_dir, sprintf('afdf%s',spm_filename))};
        matlabbatch{batch_idx}.spm.meeg.preproc.correct.mode = 'berg';
        matlabbatch{batch_idx}.spm.meeg.preproc.correct.prefix = 'T';
        batch_idx=batch_idx+1;
        spm_jobman('run',matlabbatch);
    end
    
    % Remove bad channels
    channels_removed=remove_bad_channels(subj_info, session_num, fullfile(analysis_dir, sprintf('Tafdf%s',spm_filename)));
    
    clear jobs
    matlabbatch={};
    batch_idx=1;

    if params.epoch
        
        % Epoch - aligned to instruction
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.D = {fullfile(analysis_dir, sprintf('Tafdf%s',spm_filename))};
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.timewin = params.instr_epoch;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.conditionlabel = 'instr';
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.eventtype = 'UPPT001_up';
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.eventvalue = 50;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.trlshift = 0;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.bc = 0;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.eventpadding = 0;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.prefix = fullfile(analysis_dir, 'instr_');
        batch_idx=batch_idx+1;
    
        % Epoch - aligned to response
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.D = {fullfile(analysis_dir, sprintf('Tafdf%s',spm_filename))};
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.timewin = params.resp_epoch;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.conditionlabel = 'resp';
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.eventtype = 'UPPT001_up';
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.eventvalue = 60;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.trlshift = 0;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.bc = 0;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.eventpadding = 0;
        matlabbatch{batch_idx}.spm.meeg.preproc.epoch.prefix = fullfile(analysis_dir, 'resp_');
        batch_idx=batch_idx+1;
    end
    
    matlabbatch{batch_idx}.spm.meeg.other.delete.D = {fullfile(analysis_dir, sprintf('afdf%s',spm_filename))};
    batch_idx=batch_idx+1;

    matlabbatch{batch_idx}.spm.meeg.other.delete.D = {fullfile(analysis_dir, sprintf('eyeblinkafdf%s',spm_filename))};
    batch_idx=batch_idx+1;
    
    matlabbatch{batch_idx}.spm.meeg.other.delete.D = {fullfile(analysis_dir, sprintf('fdf%s',spm_filename))};
    batch_idx=batch_idx+1;

    matlabbatch{batch_idx}.spm.meeg.other.delete.D = {fullfile(analysis_dir, sprintf('df%s',spm_filename))};
    batch_idx=batch_idx+1;

    matlabbatch{batch_idx}.spm.meeg.other.delete.D = {fullfile(analysis_dir, sprintf('f%s',spm_filename))};
    batch_idx=batch_idx+1;

    matlabbatch{batch_idx}.spm.meeg.other.delete.D = {fullfile(analysis_dir, spm_filename)};
    batch_idx=batch_idx+1;
    
    spm_jobman('run',matlabbatch);
    
    % Label trial conditions
    [no_response incorrect]=label_trials(fullfile(analysis_dir, sprintf('instr_Tafdf%s',spm_filename)), stim_file, data_file, 'delete_no_resp',true);
    [no_response incorrect]=label_trials(fullfile(analysis_dir, sprintf('resp_Tafdf%s',spm_filename)), stim_file, data_file, 'delete_no_resp',true);
    
    instr_run_files{end+1,1}=fullfile(analysis_dir, sprintf('instr_Tafdf%s',spm_filename));
    resp_run_files{end+1,1}=fullfile(analysis_dir, sprintf('resp_Tafdf%s',spm_filename));
        
    close all;

end

clear jobs
matlabbatch={};
batch_idx=1;

if subj_info.sessions(session_num)>1
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.D = instr_run_files;
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.file = '.*';
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.labelorg = '.*';
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.labelnew = '#labelorg#';
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.prefix = 'c';
    batch_idx=batch_idx+1;

    matlabbatch{batch_idx}.spm.meeg.preproc.merge.D = resp_run_files;
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.file = '.*';
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.labelorg = '.*';
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.labelnew = '#labelorg#';
    matlabbatch{batch_idx}.spm.meeg.preproc.merge.prefix = 'c';
    batch_idx=batch_idx+1;

    [path filename]=fileparts(instr_run_files{1,1});
    matlabbatch{batch_idx}.spm.meeg.other.copy.D = {sprintf('c%s.mat', filename)};
    matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = fullfile(analysis_dir, sprintf('cinstr_Tafdf%d.mat', session_num));
    batch_idx=batch_idx+1;

    matlabbatch{batch_idx}.spm.meeg.other.delete.D = {sprintf('c%s.mat', filename)};
    batch_idx=batch_idx+1;

    [path filename]=fileparts(resp_run_files{1,1});
    matlabbatch{batch_idx}.spm.meeg.other.copy.D = {sprintf('c%s.mat', filename)};
    matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = fullfile(analysis_dir, sprintf('cresp_Tafdf%d.mat', session_num));
    batch_idx=batch_idx+1;

    matlabbatch{batch_idx}.spm.meeg.other.delete.D = {sprintf('c%s.mat', filename)};
    batch_idx=batch_idx+1;
else
    matlabbatch{batch_idx}.spm.meeg.other.copy.D = {fullfile(analysis_dir, sprintf('instr_Tafdf%d-1.mat',session_num))};
    matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = fullfile(analysis_dir, sprintf('cinstr_Tafdf%d.mat', session_num));
    batch_idx=batch_idx+1;
    
    matlabbatch{batch_idx}.spm.meeg.other.copy.D = {fullfile(analysis_dir, sprintf('resp_Tafdf%d-1.mat',session_num))};
    matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = fullfile(analysis_dir, sprintf('cresp_Tafdf%d.mat', session_num));
    batch_idx=batch_idx+1;
end

spm_jobman('run',matlabbatch);


bad_trials=[];

% Artifact rejection - instruction aligned
original_data=spm_eeg_load(fullfile(analysis_dir, sprintf('cinstr_Tafdf%d.mat', session_num)));
instr_old_trials=original_data.badtrials();
instr_old_channels=original_data.badchannels();
S=[];
S.D=original_data;
S.outfile=fullfile(analysis_dir, sprintf('rcinstr_Tafdf%d.mat', session_num));
spm_eeg_copy(S);
[instr_mean_trial_var instr_std_trial_var]=plot_trial_var(original_data);
plot_channel_var(original_data);

instr=spm_eeg_load(fullfile(analysis_dir, sprintf('rcinstr_Tafdf%d.mat', session_num)));
good_chans=setdiff([33:307],instr.badchannels);
trial_var=squeeze(mean(std(instr(good_chans,:,:),[],1).^2,2));
instr_bad_trials=find(abs(trial_var-instr_mean_trial_var)>=2.5*instr_std_trial_var);
bad_trials=union(bad_trials, instr_bad_trials);

% Artifact rejection - response aligned
original_data=spm_eeg_load(fullfile(analysis_dir, sprintf('cresp_Tafdf%d.mat', session_num)));
resp_old_trials=original_data.badtrials();
resp_old_channels=original_data.badchannels();
S=[];
S.D=original_data;
S.outfile=fullfile(analysis_dir, sprintf('rcresp_Tafdf%d.mat', session_num));
spm_eeg_copy(S);
[resp_mean_trial_var resp_std_trial_var]=plot_trial_var(original_data);
plot_channel_var(original_data);

resp=spm_eeg_load(fullfile(analysis_dir, sprintf('rcresp_Tafdf%d.mat', session_num)));
good_chans=setdiff([33:307],resp.badchannels);
trial_var=squeeze(mean(std(resp(good_chans,:,:),[],1).^2,2));
resp_bad_trials=find(abs(trial_var-resp_mean_trial_var)>=2.5*resp_std_trial_var);
bad_trials=union(bad_trials,resp_bad_trials);

instr_filtered_data=badtrials(instr,bad_trials,1);
instr_filtered_data.save();
plot_trial_var(instr_filtered_data, 'mean_trial_var', instr_mean_trial_var, 'std_trial_var', instr_std_trial_var);
plot_channel_var(instr_filtered_data);
artifact_trials=setdiff(instr_filtered_data.badtrials(),instr_old_trials);
artifact_channels=setdiff(instr_filtered_data.badchannels(),instr_old_channels);
removed_trials={};
for idx=1:length(artifact_trials)
    removed_trials{end+1}=num2str(artifact_trials(idx));
end
removed_channels={};
for idx=1:length(artifact_channels)
    removed_channels{end+1}=S.D.chanlabels{artifact_channels(idx)};
end


resp_filtered_data=badtrials(resp,bad_trials,1);
resp_filtered_data.save();
plot_trial_var(resp_filtered_data, 'mean_trial_var', resp_mean_trial_var, 'std_trial_var', resp_std_trial_var);
plot_channel_var(resp_filtered_data);
artifact_trials=setdiff(resp_filtered_data.badtrials(),resp_old_trials);
artifact_channels=setdiff(resp_filtered_data.badchannels(),resp_old_channels);
removed_trials={};
for idx=1:length(artifact_trials)
    removed_trials{end+1}=num2str(artifact_trials(idx));
end
removed_channels={};
for idx=1:length(artifact_channels)
    removed_channels{end+1}=S.D.chanlabels{artifact_channels(idx)};
end


clear jobs
matlabbatch={};
batch_idx=1;

% Remove prefiltered files
matlabbatch{batch_idx}.spm.meeg.other.delete.D = {fullfile(analysis_dir, sprintf('cresp_Tafdf%d.mat', session_num))};
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.meeg.other.delete.D = {fullfile(analysis_dir, sprintf('cinstr_Tafdf%d.mat', session_num))};
batch_idx=batch_idx+1;

% Remove run epoched files1
matlabbatch{batch_idx}.spm.meeg.other.delete.D = instr_run_files;
batch_idx=batch_idx+1;
    
matlabbatch{batch_idx}.spm.meeg.other.delete.D = resp_run_files;
batch_idx=batch_idx+1;

% Average
matlabbatch{batch_idx}.spm.meeg.averaging.average.D = {fullfile(analysis_dir, sprintf('rcinstr_Tafdf%d.mat', session_num))};
matlabbatch{batch_idx}.spm.meeg.averaging.average.userobust.standard = false;
matlabbatch{batch_idx}.spm.meeg.averaging.average.plv = false;
matlabbatch{batch_idx}.spm.meeg.averaging.average.prefix = fullfile(analysis_dir, 'm');
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.meeg.averaging.average.D = {fullfile(analysis_dir, sprintf('rcresp_Tafdf%d.mat', session_num))};
matlabbatch{batch_idx}.spm.meeg.averaging.average.userobust.standard = false;
matlabbatch{batch_idx}.spm.meeg.averaging.average.plv = false;
matlabbatch{batch_idx}.spm.meeg.averaging.average.prefix = fullfile(analysis_dir, 'm');
batch_idx=batch_idx+1;

spm_jobman('run',matlabbatch);

close all;
