function plot_channel_var(D, varargin)

defaults = struct('output_file','','output_format','png');  %define default values
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
chan_var=squeeze(mean(std(D(good_chans,:,good_trials),[],3).^2,2));
mean_chan_var=mean(chan_var);
std_chan_var=std(chan_var);
hold all;
plot(good_chans,chan_var,'o')
plot([min(good_chans) max(good_chans)], [mean_chan_var mean_chan_var],'k--');
%plot([min(good_chans) max(good_chans)], [mean_chan_var-2.5*std_chan_var mean_chan_var-2.5*std_chan_var],'r');
%plot([min(good_chans) max(good_chans)], [mean_chan_var+2.5*std_chan_var mean_chan_var+2.5*std_chan_var],'r');
xlabel('Channel');
ylabel('Variance');

if length(params.output_file)>0
    saveas(fig, params.output_file, params.output_format);
    close(fig);
end
