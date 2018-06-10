function subjects_behavior(subjects, data_dir, varargin)

defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

congruence_conditions={'congruent','incongruent'};
coherence_conditions={'low','med','high'};
conditions={'congruent - low','congruent - med','congruent - high',...
    'incongruent - low','incongruent - med','incongruent - high'};
rts=dict();
correct=dict();

for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    [subject_rts, subject_correct]=subject_behavior(subj_info, fullfile(data_dir, subj_info.subj_id), 'plot',false);
    for cond_idx=1:length(conditions)
        condition=conditions{cond_idx};
        
        x=rts(condition);
        x(end+1)=mean(subject_rts(condition));
        rts(condition)=x;
        
        x=correct(condition);
        x(end+1)=mean(subject_correct(condition));
        correct(condition)=x;        
    end
end

mean_rts=[mean(rts('congruent - low')) mean(rts('incongruent - low'));...
    mean(rts('congruent - med')) mean(rts('incongruent - med'));...
    mean(rts('congruent - high')) mean(rts('incongruent - high'))];
std_rts=[std(rts('congruent - low'))/sqrt(length(rts('congruent - low'))) std(rts('incongruent - low'))/sqrt(length(rts('incongruent - low')));...
    std(rts('congruent - med'))/sqrt(length(rts('congruent - med'))) std(rts('incongruent - med'))/sqrt(length(rts('incongruent - med')));...
    std(rts('congruent - high'))/sqrt(length(rts('congruent - high'))) std(rts('incongruent - low'))/sqrt(length(rts('incongruent - high')))];

mean_correct=[mean(correct('congruent - low')) mean(correct('incongruent - low'));...
    mean(correct('congruent - med')) mean(correct('incongruent - med'));...
    mean(correct('congruent - high')) mean(correct('incongruent - high'))];
std_correct=[std(correct('congruent - low'))/sqrt(length(correct('congruent - low'))) std(correct('incongruent - low'))/sqrt(length(correct('incongruent - low')));...
    std(correct('congruent - med'))/sqrt(length(correct('congruent - med'))) std(correct('incongruent - med'))/sqrt(length(correct('incongruent - med')));...
    std(correct('congruent - high'))/sqrt(length(correct('congruent - high'))) std(correct('incongruent - low'))/sqrt(length(correct('incongruent - high')))];

fig=figure();
h=barwitherr(std_rts, mean_rts);
set(gca, 'XTickLabel', {'low', 'medium', 'high'});
legend('congruent','incongruent');
xlabel('Coherence');
ylabel('RT (ms)');

fig=figure();
h=barwitherr(std_correct.*100, mean_correct.*100);
set(gca, 'XTickLabel', {'low', 'medium', 'high'});
legend('congruent','incongruent');
xlabel('Coherence');
ylabel('% correct');
