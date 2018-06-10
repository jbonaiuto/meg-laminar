function [rts, correct]=run_behavior(subj_info, session_num, run_num, data_dir, varargin)

defaults = struct('plot',true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

load(fullfile(data_dir, 'behavior', sprintf('data_%d.mat', run_num)));

load(fullfile(data_dir, 'behavior', sprintf('stim_%d.mat', run_num)));

% Correct for left/right mismatch
stim.trials(:,1)=1+2-stim.trials(:,1);
stim.trials(:,4)=1+2-stim.trials(:,4);

low_congruence=.5*stim.threshold;
med_congruence=stim.threshold;
high_congruence=1.5*stim.threshold;

correct_trials=data.responses(:,1)==stim.trials(:,4);

cong_low_trials=find(stim.trials(:,3)==1 & stim.trials(:,2)==low_congruence);
cong_med_trials=find(stim.trials(:,3)==1 & stim.trials(:,2)==med_congruence);
cong_high_trials=find(stim.trials(:,3)==1 & stim.trials(:,2)==high_congruence);

incong_low_trials=find(stim.trials(:,3)==0 & stim.trials(:,2)==low_congruence);
incong_med_trials=find(stim.trials(:,3)==0 & stim.trials(:,2)==med_congruence);
incong_high_trials=find(stim.trials(:,3)==0 & stim.trials(:,2)==high_congruence);

rts=dict();
rts('congruent - low')=data.responses(cong_low_trials,2).*1000;
rts('congruent - med')=data.responses(cong_med_trials,2).*1000;
rts('congruent - high')=data.responses(cong_high_trials,2).*1000;
rts('incongruent - low')=data.responses(incong_low_trials,2).*1000;
rts('incongruent - med')=data.responses(incong_med_trials,2).*1000;
rts('incongruent - high')=data.responses(incong_high_trials,2).*1000;

correct=dict();
correct('congruent - low')=correct_trials(cong_low_trials);
correct('congruent - med')=correct_trials(cong_med_trials);
correct('congruent - high')=correct_trials(cong_high_trials);
correct('incongruent - low')=correct_trials(incong_low_trials);
correct('incongruent - med')=correct_trials(incong_med_trials);
correct('incongruent - high')=correct_trials(incong_high_trials);

if params.plot
    mean_rts=[nanmean(rts('congruent - low')) nanmean(rts('incongruent - low'));...
        nanmean(rts('congruent - med')) nanmean(rts('incongruent - med'));...
        nanmean(rts('congruent - high')) nanmean(rts('incongruent - high'))];
    std_rts=[nanstd(rts('congruent - low'))/sqrt(length(rts('congruent - low'))) nanstd(rts('incongruent - low'))/sqrt(length(rts('incongruent - low')));...
        nanstd(rts('congruent - med'))/sqrt(length(rts('congruent - med'))) nanstd(rts('incongruent - med'))/sqrt(length(rts('incongruent - med')));...
        nanstd(rts('congruent - high'))/sqrt(length(rts('congruent - high'))) nanstd(rts('incongruent - low'))/sqrt(length(rts('incongruent - high')))];

    mean_correct=[mean(correct('congruent - low')) mean(correct('incongruent - low'));...
        mean(correct('congruent - med')) mean(correct('incongruent - med'));...
        mean(correct('congruent - high')) mean(correct('incongruent - high'))];
    std_correct=[std(correct('congruent - low'))/sqrt(length(correct('congruent - low'))) std(correct('incongruent - low'))/sqrt(length(correct('incongruent - low')));...
        std(correct('congruent - med'))/sqrt(length(correct('congruent - med'))) std(correct('incongruent - med'))/sqrt(length(correct('incongruent - med')));...
        std(correct('congruent - high'))/sqrt(length(correct('congruent - high'))) std(correct('incongruent - low'))/sqrt(length(correct('incongruent - high')))];

    fig=figure();
    barwitherr(std_rts, mean_rts);
    set(gca, 'XTickLabel', {'low', 'medium', 'high'});
    legend('congruent','incongruent');
    xlabel('Coherence');
    ylabel('RT (ms)');

    fig=figure();
    barwitherr(std_correct.*100, mean_correct.*100);
    set(gca, 'XTickLabel', {'low', 'medium', 'high'});
    legend('congruent','incongruent');
    xlabel('Coherence');
    ylabel('% correct');
end
