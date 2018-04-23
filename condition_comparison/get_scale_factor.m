function scale_factor=get_scale_factor(conditions, cond_trials, contrast)

cond_vals=zeros(1,length(conditions));
for cond_idx=1:length(conditions)
    condition=conditions{cond_idx};
    cond_trial=mean(cond_trials(condition),2);
    switch contrast.direction
        case 'negative'
            cond_vals(cond_idx)=min(cond_trial);
        case 'positive'
            cond_vals(cond_idx)=max(cond_trial);
    end
end
scale_factor=max(abs(cond_vals));
