function [mean_trial_var std_trial_var]=plot_trial_var(D, varargin)

defaults = struct('output_file','','output_format','png','mean_trial_var',[],'std_trial_var',[]);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

ntrials=size(D,3);
good_trials=setdiff([1:ntrials],D.badtrials);
meg_chans=[33:307];
good_chans=setdiff(meg_chans,D.badchannels);
fig=figure();
trial_var=squeeze(mean(std(D(good_chans,:,good_trials),[],1).^2,2));
if length(params.mean_trial_var)
    mean_trial_var=params.mean_trial_var;
    std_trial_var=params.std_trial_var;
else
    mean_trial_var=mean(trial_var);
    std_trial_var=std(trial_var);
end
hold all
plot(good_trials,trial_var,'o');
plot([min(good_trials) max(good_trials)],[mean_trial_var mean_trial_var],'k--');
plot([min(good_trials) max(good_trials)],[mean_trial_var-2.5*std_trial_var mean_trial_var-2.5*std_trial_var],'r');
plot([min(good_trials) max(good_trials)],[mean_trial_var+2.5*std_trial_var mean_trial_var+2.5*std_trial_var],'r');
hold off;
xlabel('Trial');
ylabel('Variance');

if length(params.output_file)>0
    saveas(fig, params.output_file, params.output_format);
    close(fig);
end
