function scale_factor=get_scale_factor_woi(conditions, cond_trials_woi, contrast)

cond_vals=zeros(1,length(conditions));
for cond_idx=1:length(conditions)
    condition=conditions{cond_idx};
    cond_vals(cond_idx)=mean(cond_trials_woi(condition));    
end
switch contrast.direction
    case 'negative'
        scale_factor=abs(min(cond_vals));
    case 'positive'
        scale_factor=max(cond_vals);
end
