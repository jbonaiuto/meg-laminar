function fig=plot_power_woi(trials, conditions, labels, varargin)
    % Parse inputs
    defaults = struct('time_limits', [], 'ax', 0);  %define default values
    params = struct(varargin{:});
    for f = fieldnames(defaults)',
        if ~isfield(params, f{1}),
            params.(f{1}) = defaults.(f{1});
        end
    end
    
    fig=0;
    if params.ax==0
        fig=figure();
    end
    bar_color=[204 204 204]./255;
    
    stderrs=zeros(1,length(conditions));
    cond_y=zeros(1,length(conditions));
    for i=1:length(conditions)
        condition=conditions{i};
        n_trials=size(trials(condition),2);
        cond_y(i)=mean(trials(condition));
        stderrs(i)=std(trials(condition))./sqrt(n_trials);            
    end   
    [hBar hErrorbar]=barwitherr(stderrs,cond_y);
    set(hBar,'FaceColor', bar_color);
    set(get(get(hBar,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(hErrorbar,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    bar_width=get(hBar,'BarWidth');
    hold all
    vals=[];
    if length(labels)>0        
        for t_idx=1:length(trials(conditions{1}))
            cond_x=[1:length(conditions)];
            cond_y=[1:length(conditions)];
            for i=1:length(conditions)
                condition=conditions{i};
                cond_trials=trials(condition);
                cond_y(i)=cond_trials(t_idx);
                vals(end+1)=cond_y(i);
            end
            plot(cond_x, cond_y, 'o--');            
        end
        val_range=max(vals)-min(vals);
        ylim([min(vals)-.1*val_range max(vals)+.1*val_range]);
    else
        val_range=max(cond_y+stderrs)-min(cond_y-stderrs);
        ylim([min(cond_y-stderrs)-.1*val_range max(cond_y+stderrs)+.1*val_range]);
    end
    
    if length(labels)>0
        legend(labels,'Location','eastoutside');
    end
    hold off;
    set(gca, 'XTick', [1:length(conditions)]);
    set(gca, 'XTickLabel', conditions);
    xlim([1-bar_width/2-(1-bar_width) length(conditions)+bar_width/2+(1-bar_width)]);
    ylabel('Power');
end