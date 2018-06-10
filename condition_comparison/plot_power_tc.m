function fig=plot_power_tc(times, trials, conditions, varargin)
    % Parse inputs
    defaults = struct('time_limits', [], 'ax', 0);  %define default values
    params = struct(varargin{:});
    for f = fieldnames(defaults)',
        if ~isfield(params, f{1}),
            params.(f{1}) = defaults.(f{1});
        end
    end
    
    time_idx=[1:length(times)];
    if length(params.time_limits)>0
        time_idx=intersect(find(times>=params.time_limits(1)),find(times<=params.time_limits(2)));
    end
    
    colors={'b','r','g','m','c','k'};
    fig=0;
    if params.ax==0
        fig=figure();
    end
    hold all;
    labels={};
    for i=1:length(conditions)
        condition=conditions{i};
        n_trials=size(trials(condition),2);
        if n_trials>1
            mean_trials=mean(trials(condition),2);
            stderr_trials=std(trials(condition),0,2)/sqrt(n_trials);
            H=shadedErrorBar(times(time_idx),mean_trials(time_idx),stderr_trials(time_idx),sprintf('%s',colors{i}));
            set(get(get(H.patch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            set(get(get(H.edge(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            set(get(get(H.edge(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            set(get(get(H.mainLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
            labels{end+1}=condition;
        else
            mean_trials=trials(condition);
            plot(times(time_idx), mean_trials(time_idx), sprintf('%s',colors{i}));
            labels{end+1}=condition;
        end
    end
    hold off;
    if length(params.time_limits)>0
        xlim(params.time_limits);
    else
        xlim([times(1) times(end)]);
    end
    legend(labels,'location','northwest');
    ylabel('Power');
    xlabel('time (ms)');
end