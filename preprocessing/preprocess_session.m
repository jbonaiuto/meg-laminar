function preprocess_session(subj_info, session_num, varargin)

defaults = struct('data_dir','/data/pred_coding','delete_last_resp',false,'highpass_freq',2.0,'downsample',250,'lowpass_freq',100.0,'instr_epoch',[-3500 1500],'resp_epoch',[-2000 2000], ...
    'blink_channel','MLT31','manual_reject',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Add path to report generation code
addpath /home/jbonaiuto/Documents/MATLAB/m2html/
    
% Directory containing session data
session_dir=fullfile(params.data_dir,'scanning', subj_info.subj_id, num2str(session_num));
% Directory to put results
analysis_dir=fullfile(params.data_dir,'analysis',subj_info.subj_id,num2str(session_num));
if exist(analysis_dir,'dir')~=7
    mkdir(analysis_dir);
end
% Directory to put report
report_dir=fullfile(analysis_dir,'preprocessing');
if exist(report_dir,'dir')~=7
    mkdir(report_dir);
end

spm('defaults','eeg');
spm_jobman('initcfg');

run_files={};
for run_num=1:subj_info.sessions(session_num)
    
    run_tpl = template('/home/jbonaiuto/Projects/meg/pred_coding/src/matlab/analysis/preprocessing/templates/','keep');
    run_tpl = set(run_tpl,'file',{'run'},{'run_template.tpl'});
    run_tpl = set(run_tpl,'var',{'PAGETITLE'},{sprintf('%s%s: Session %d - Run %d', subj_info.subj_id, subj_info.birth_date, session_num, run_num)});


    % Directory containing run data
    run_code=[subj_info.subj_id subj_info.birth_date '_JamesBonaiuto_' subj_info.scan_date{session_num} '_0' num2str(run_num)];
    run_dir=fullfile(session_dir, [run_code '.ds']);

    % File containing stimulus information for each trial
    stim_file=fullfile(session_dir, ['stim_' subj_info.subj_id '_' num2str(run_num) '.mat']);
    % File containing behavioral data from each trial
    data_file=fullfile(session_dir, ['data_' subj_info.subj_id '_' num2str(run_num) '.mat']);

    spm_filename=sprintf('%s-%d-%d.mat',subj_info.subj_id,session_num,run_num);
    
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
    
    file_path=fullfile('analysis',subj_info.subj_id,num2str(session_num),'preprocessing',sprintf('fiducial_coils-%d.png',run_num));
    plot_fiducial_coils(fullfile(analysis_dir, spm_filename), 'lim_last_event',false,'lim_jump',false,'output_file', fullfile(params.data_dir,file_path));
    run_tpl = set(run_tpl, 'var', {'FIDUCIALSRC'}, {['http://fortressofjollitude.zapto.org/' file_path]});

    file_path=fullfile('analysis',subj_info.subj_id,num2str(session_num),'preprocessing',sprintf('fiducial_coils_lastevent-%d.png',run_num));
    plot_fiducial_coils(fullfile(analysis_dir, spm_filename), 'lim_last_event',true,'lim_jump',false, 'output_file', fullfile(params.data_dir,file_path));
    run_tpl = set(run_tpl, 'var', {'FIDUCIALLASTEVENTSRC'}, {['http://fortressofjollitude.zapto.org/' file_path]});

    file_path=fullfile('analysis',subj_info.subj_id,num2str(session_num),'preprocessing',sprintf('fiducial_coils_nojump-%d.png',run_num));
    plot_fiducial_coils(fullfile(analysis_dir, spm_filename), 'lim_last_event', false, 'lim_jump', true, 'output_file', fullfile(params.data_dir,file_path));
    run_tpl = set(run_tpl, 'var', {'FIDUCIALNOJUMPSRC'}, {['http://fortressofjollitude.zapto.org/' file_path]});

    % Correct event timings, adjust initial events
    file_path=fullfile('analysis',subj_info.subj_id,num2str(session_num),'preprocessing',sprintf('diode-%d.png',run_num));
    [num_diode_onsets num_dots_evts num_instr_evts num_resp_evts]=adjust_run_events(subj_info, session_num, run_num, 'data_dir', params.data_dir', 'plot_diode',true,...
        'output_file', fullfile(params.data_dir,file_path), 'delete_last_resp',params.delete_last_resp,'delete_no_resp',true);
    run_tpl = set(run_tpl,'var',{'DIODECHANNEL','DIODEONSETS','DIODETHRESH','DIODESRC','DOTSEVENTS','INSTREVENTS','RESPEVENTS'},{num2str(subj_info.diode_ch(session_num)),
        num2str(num_diode_onsets),num2str(subj_info.diode_thresh(session_num)),['http://fortressofjollitude.zapto.org/' file_path],
        num2str(num_dots_evts),num2str(num_instr_evts),num2str(num_resp_evts)});

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
    run_tpl = set(run_tpl,'var',{'HIGHPASSFREQ'},{num2str(params.highpass_freq)});

    % Downsample
    matlabbatch{batch_idx}.spm.meeg.preproc.downsample.D = {fullfile(analysis_dir, sprintf('f%s',spm_filename))};
    matlabbatch{batch_idx}.spm.meeg.preproc.downsample.fsample_new = params.downsample;
    matlabbatch{batch_idx}.spm.meeg.preproc.downsample.prefix = fullfile(analysis_dir,'d');
    batch_idx=batch_idx+1;
    run_tpl = set(run_tpl,'var',{'DOWNSAMPLE'},{num2str(params.downsample)});

    % Lowpass filter
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.D = {fullfile(analysis_dir, sprintf('df%s',spm_filename))};
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.type = 'butterworth';
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.band = 'low';
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.freq = params.lowpass_freq;
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.dir = 'twopass';
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.order = 5;
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.prefix = fullfile(analysis_dir,'f');
    batch_idx=batch_idx+1;
    run_tpl = set(run_tpl,'var',{'LOWPASSFREQ'},{num2str(params.lowpass_freq)});
    spm_jobman('run',matlabbatch);
    
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

    matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.D = {fullfile(analysis_dir, sprintf('eyeblinkafdf%s',spm_filename))};
    matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.mode{1}.svd.timewin = [-Inf Inf];
    matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.mode{1}.svd.threshold = NaN;
    matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.mode{1}.svd.ncomp = 4;
    batch_idx=batch_idx+1;
    
    spm_jobman('run',matlabbatch);
    
    clear jobs
    matlabbatch={};
    batch_idx=1;
    
    ncomp=prompt('Number of components: ');
    matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.D = {fullfile(analysis_dir, sprintf('eyeblinkafdf%s',spm_filename))};
    matlabbatch{batch_idx}.spm.meeg.preproc.sconfounds.mode{1}.clear=1;
    
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
    
    run_tpl = set(run_tpl,'var',{'BLINKSREMOVED'},{['TRUE (' params.blink_channel ')']});
    Fgraph = spm_figure('GetWin','Graphics');
    file_path=fullfile('analysis',subj_info.subj_id,num2str(session_num),'preprocessing',sprintf('blink_component-%d.png',run_num));
    saveas(Fgraph, fullfile(params.data_dir, file_path));
    run_tpl = set(run_tpl,'var',{'BLINKCMPNTSRC'},{['http://fortressofjollitude.zapto.org/' file_path]});
    
    % Remove bad channels
    channels_removed=remove_bad_channels(subj_info, session_num, fullfile(analysis_dir, sprintf('Tafdf%s',spm_filename)));
    run_tpl = set(run_tpl,'var',{'CHANNELSREMOVED'},{strjoin(channels_removed, ', ')});
    
    clear jobs
    matlabbatch={};
    batch_idx=1;

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
    run_tpl = set(run_tpl,'var',{'INSTREPOCH'},{[num2str(params.instr_epoch(1)) '-' num2str(params.instr_epoch(2)) 'ms']});

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
    run_tpl = set(run_tpl,'var',{'RESPEPOCH'},{[num2str(params.resp_epoch(1)) '-' num2str(params.resp_epoch(2)) 'ms']});
    
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
    run_tpl = set(run_tpl,'var',{'NORESP','INCORRECT'},{num2str(no_response),num2str(incorrect)});

    instr_run_files{end+1,1}=sprintf('instr_Tafdf%s',spm_filename);
    resp_run_files{end+1,1}=sprintf('resp_Tafdf%s',spm_filename);
    
    run_tpl = parse(run_tpl,'OUT', {'run'});

    %- Display the content of the output of the parsing
    out_file=fullfile(report_dir,sprintf('preprocessing-%d.html',run_num));
    fid=fopen(out_file,'w');
    %display(get(tpl,'OUT'))
    fprintf(fid, strrep(get(run_tpl,'OUT'),'%','%%'));
    fclose(fid);

    close all;

end

clear jobs
matlabbatch={};
batch_idx=1;

matlabbatch{batch_idx}.spm.meeg.preproc.merge.D = instr_run_files;
matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.file = '.*';
matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.labelorg = '.*';
matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.labelnew = '#labelorg#';
matlabbatch{batch_idx}.spm.meeg.preproc.merge.prefix = fullfile(analysis_dir,'c');
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.meeg.preproc.merge.D = resp_run_files;
matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.file = '.*';
matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.labelorg = '.*';
matlabbatch{batch_idx}.spm.meeg.preproc.merge.recode.labelnew = '#labelorg#';
matlabbatch{batch_idx}.spm.meeg.preproc.merge.prefix = fullfile(analysis_dir,'c');
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.meeg.other.copy.D = {fullfile(analysis_dir, sprintf('c%s', instr_run_files{1,1}))};
matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = fullfile(analysis_dir, sprintf('cinstr_Tafdf%d-%d.mat', subj_info.subj_id, session_num));
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.meeg.other.delete.D = {fullfile(analysis_dir, sprintf('c%s', instr_run_files{1,1}))};
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.meeg.other.copy.D = {fullfile(analysis_dir, sprintf('c%s', resp_run_files{1,1}))};
matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = fullfile(analysis_dir, sprintf('cresp_Tafdf%d-%d.mat', subj_info.subj_id, session_num));
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.meeg.other.delete.D = {fullfile(analysis_dir, sprintf('c%s', resp_run_files{1,1}))};
batch_idx=batch_idx+1;

spm_jobman('run',matlabbatch);