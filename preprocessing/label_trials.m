function [no_response incorrect]=label_trials(meg_file, stim_file, data_file, varargin)
% function label_trials(meg_file, stim_file, data_file, varargin)
% Label trial conditions
% INPUT: 
%   meg_file: file containing MEG data
%   stim_file: file containing stimulus information for each trial
%   data_file: file containing behavioral data for each trial
%   response_aligned: whether or not epochs are aligned to response (default=false)
%   label_congruence: whether or not to label trial congruence (default=true)
% ---------------------------
% JJB (j.bonaiuto@ucl.ac.uk) Jul 2016
% 

% Parse inputs
defaults = struct('delete_no_resp',false,'label_congruence',true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end
    
% Load files
load(meg_file);
load(stim_file);
load(data_file);

% Correct stim direction mismatch
stim.trials(:,1)=1+2-stim.trials(:,1);
stim.trials(:,4)=1+2-stim.trials(:,4);
% Compute trial accuracy
data.correct=stim.trials(:,4)==data.responses(:,1);

% Compute coherence levels
low_coherence=min(stim.trials(:,2));
med_coherence=low_coherence*2.0;
high_coherence=low_coherence*3.0;

% Construct list of conditions
if params.label_congruence
    D.condlist={'congruent - low','congruent - med','congruent - high','incongruent - low','incongruent - med','incongruent - high'};
else
    D.condlist={'low','med','high'};
end

no_response=0;
incorrect=0;
for behav_idx=1:length(stim.trials)
    meg_idx=behav_idx;
    % Reduce MEG idx by number of no responses (these trials will not be in MEG dataset after epoching to response event)
    if params.delete_no_resp
        meg_idx=meg_idx-no_response;
    end

    cond_label='';
    % Prepend condition label with incongruent or congruent if labeling congruence
    % Append coherence level to end of condition label
    if stim.trials(behav_idx,2)==low_coherence
        cond_label='low';
    elseif stim.trials(behav_idx,2)==med_coherence
        cond_label='med';
    elseif stim.trials(behav_idx,2)==high_coherence
        cond_label='high';
    else
        ['Error! coherence=' num2str(stim.trials(behav_idx,2))]
    end        
    if params.label_congruence
        if stim.trials(behav_idx,3)==0
            cond_label=['incongruent - ' cond_label];
        else
            cond_label=['congruent - ' cond_label];
        end
    end
    if meg_idx<=length(D.trials)
        D.trials(meg_idx).label=cond_label;
    end
    % No response
    if data.responses(behav_idx,1)==0
        no_response=no_response+1;
    % Incorrect
    elseif data.correct(behav_idx)==0
        incorrect=incorrect+1;
    end
end

disp(['no_response=' num2str(no_response) ', incorrect=' num2str(incorrect)]); 
save(meg_file,'D');

