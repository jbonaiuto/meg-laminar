function [rts, correct]=session_behavior(subj_info, session_num, varargin)

defaults = struct('plot',true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

conditions={'congruent - low','congruent - med','congruent - high',...
    'incongruent - low','incongruent - med','incongruent - high'};
rts=dict();
correct=dict();

for run_num=1:subj_info.sessions(session_num)
    [run_rts, run_correct]=run_behavior(subj_info, session_num, run_num,...
        'plot',false);
    for cond_idx=1:length(conditions)
        condition=conditions{cond_idx};
        
        x=rts(condition);
        x(end+1:end+length(run_rts(condition)))=run_rts(condition);
        rts(condition)=x;
        
        x=correct(condition);
        x(end+1:end+length(run_correct(condition)))=run_correct(condition);
        correct(condition)=x;        
    end
end

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