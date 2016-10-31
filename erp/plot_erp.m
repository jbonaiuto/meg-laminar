function plot_erp(trial_times, trial_erps, conditions, varargin)

defaults = struct('time_limits',[],'output_file','','output_format','png');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

colors={'b','r','g','b--','r--','g--'};
fig=figure('Position',[1 1 800 600],'PaperUnits','points','PaperPosition',[1 1 800 600],'PaperPositionMode','manual');
hold all;
for i=1:length(conditions)
    condition=conditions{i};
    n_trials=size(trial_erps(condition),2);
    if n_trials>1
        mean_erp=mean(trial_erps(condition),2);
        stderr_erp=std(trial_erps(condition),0,2)/sqrt(n_trials);
        H=shadedErrorBar(trial_times,mean_erp,stderr_erp,colors{i});
        set(get(get(H.patch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(get(get(H.edge(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(get(get(H.edge(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(get(get(H.mainLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
    else
        plot(trial_times, trial_erps(condition), colors{i});
    end
end
hold off;
if length(params.time_limits)>0
    xlim(params.time_limits);
end
legend(conditions);
ylabel('field intensity (fT)');
xlabel('time (ms)');
%title(['channel ' chan_label ' (MEGGRAD)']);
if length(params.output_file)>0
    saveas(fig, params.output_file, params.output_format);
    close(fig);
end
