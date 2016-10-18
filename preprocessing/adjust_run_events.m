function [num_diode_onsets num_dots_evts num_instr_evts num_resp_evts]=adjust_run_events(subj_info, session_num, run_num, varargin)
% function adjust_run_events(subj_info, session_num, run_num, varargin)
% Adjust timing of events in a run
% INPUT: 
%   subj_info: subject info structure
%   session_num: session number
%   run_num: run number
%   diode_ch: channel containing the diode signal (default=316)
% ---------------------------
% JJB (j.bonaiuto@ucl.ac.uk) Jul 2016
% 

% Parse inputs
defaults = struct('data_dir','/data/pred_coding','plot_diode',false,'output_file','','output_format','png','delete_last_resp',false,'delete_no_resp',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% File containing data
analysis_dir=fullfile(params.data_dir,'analysis',subj_info.subj_id,num2str(session_num));
spm_file_name=fullfile(analysis_dir, sprintf('%s-%d-%d.mat',subj_info.subj_id,session_num,run_num));

% Adjust event timings
load(spm_file_name);
% Starts four time steps after last initial trigger
start_time=find(D.data(1,:)==70,1)+4;
% End time
end_time=size(D.data,2);

% Use subject-specific threshold to find when diode in in up state
diode_up_times=find(D.data(subj_info.diode_ch(session_num),:)>subj_info.diode_thresh(session_num));
% Only use up state times from after the start of the run
diode_up_times=diode_up_times(find(diode_up_times>start_time));
% Find changes in diode state
diode_diff=diode_up_times(2:length(diode_up_times))-diode_up_times(1:length(diode_up_times)-1);
diode_times=[];
if length(diode_up_times)>0
    % Find where diode goes from down to up
    diode_times=[diode_up_times(1) diode_up_times(find(diode_diff>1)+1)];
    if length(diode_times)==361 || length(diode_times)==181
        diode_times=diode_times(1:end-1);
    end
end
if params.plot_diode
    fig=figure('Position',[1 1 1800 400],'PaperUnits','points','PaperPosition',[1 1 900 200],'PaperPositionMode','manual');
    subplot(1,1,1);
    hold on
    plot(D.data(subj_info.diode_ch(session_num),:),'b');
    plot([1 size(D.data,2)],[subj_info.diode_thresh(session_num) subj_info.diode_thresh(session_num)],'r');
    for i=1:length(diode_times)
        plot([diode_times(i) diode_times(i)],[min(D.data(subj_info.diode_ch(session_num),:)) subj_info.diode_thresh(session_num)],'g');
    end
    hold off
    if length(params.output_file)>0
        saveas(fig, params.output_file, params.output_format);
    end
end
num_diode_onsets=length(diode_times);
disp(['num diode onsets=' num2str(num_diode_onsets)]);

% In early sessions there was only a diode signal for the dots
dots_onset=diode_times;
% In later sessions there was a diode for the dots and the instruction stimulus
if length(diode_times)>180
    dots_onset=diode_times(1:2:end);
    instr_onset=diode_times(2:2:end);
end

% Count events
first_dots_idx=0;
second_dots_idx=0;
num_dots_evts=0;
for i=1:length(D.trials.events)
    % dots event
    if strcmp(D.trials.events(i).type,'UPPT001_up') && D.trials.events(i).value==30
        if first_dots_idx==0
            first_dots_idx=i;
        elseif second_dots_idx==0
            second_dots_idx=i;
        end
        num_dots_evts=num_dots_evts+1;
    end
end
% Remove first events if more than 180 of them
if num_dots_evts>180
    % Delete all events up to second trial event
    evts_to_delete=[1:second_dots_idx-1];
    D.trials.events=D.trials.events(setdiff([1:length(D.trials.events)],evts_to_delete));        
end

if params.delete_last_resp
    % Delete last response
    last_resp=-1;
    for i=1:length(D.trials.events)
        if strcmp(D.trials.events(i).type,'UPPT001_up') && D.trials.events(i).value==60
            last_resp=i;
        end
    end
    D.trials.events=D.trials.events(setdiff([1:length(D.trials.events)],[last_resp]));        
end

trial_idx=1;
for i=1:length(D.trials.events)
    % If this is a trigger event that occurs after the run starts
    if strcmp(D.trials.events(i).type,'UPPT001_up')
        % Correct dots onset
        if D.trials.events(i).value==30
            % Use photodiode signal if exists
            if length(diode_times)==180 || length(diode_times)==360
                D.trials.events(i).time=dots_onset(trial_idx)/D.Fsample;
            % Use mean delay if not
            else
                D.trials.events(i).time=D.trials.events(i).time+0.0301856;
            end
        % Correct instruction onset
        elseif D.trials.events(i).value==50
            % Use photodiode signal if exists
            if length(diode_times)==360
                D.trials.events(i).time=instr_onset(trial_idx)/D.Fsample;
            % Use mean delay if not
            else
                % Mean delay is different for first trial
                if trial_idx==1
                    D.trials.events(i).time=D.trials.events(i).time+0.0190278;
                else
                    D.trials.events(i).time=D.trials.events(i).time+0.0302645;
                end
            end
            trial_idx=trial_idx+1;            
        end
    end
end

if params.delete_no_resp
    load(fullfile('/data','pred_coding','scanning', subj_info.subj_id, num2str(session_num), ['data_' subj_info.subj_id '_' num2str(run_num) '.mat']));
    evts_to_delete=[];
    no_resp_trials=find(data.responses(:,1)==0);
    num_dots_evts=0;
    num_instr_evts=0;
    for i=1:length(D.trials.events)
        % Dots event
        if strcmp(D.trials.events(i).type,'UPPT001_up') && D.trials.events(i).value==30
            num_dots_evts=num_dots_evts+1;
            if length(find(no_resp_trials==num_dots_evts))
                evts_to_delete=[evts_to_delete i];
            end
        % Instr event
        elseif strcmp(D.trials.events(i).type,'UPPT001_up') && D.trials.events(i).value==50
            num_instr_evts=num_instr_evts+1;
            if length(find(no_resp_trials==num_instr_evts))
                evts_to_delete=[evts_to_delete i];
            end
        end
    end
    D.trials.events=D.trials.events(setdiff([1:length(D.trials.events)],evts_to_delete));
end
    
num_dots_evts=0;
num_instr_evts=0;
num_resp_evts=0;
for i=1:length(D.trials.events)
    % Dots event
    if strcmp(D.trials.events(i).type,'UPPT001_up') && D.trials.events(i).value==30
        num_dots_evts=num_dots_evts+1;
    % Instr event
    elseif strcmp(D.trials.events(i).type,'UPPT001_up') && D.trials.events(i).value==50
        num_instr_evts=num_instr_evts+1;
    % Resp event
    elseif strcmp(D.trials.events(i).type,'UPPT001_up') && D.trials.events(i).value==60
        num_resp_evts=num_resp_evts+1;
    end
end
disp(['num dots events=' num2str(num_dots_evts)]);
disp(['num instr events=' num2str(num_instr_evts)]);
disp(['num resp events=' num2str(num_resp_evts)]);

save(spm_file_name,'D');
