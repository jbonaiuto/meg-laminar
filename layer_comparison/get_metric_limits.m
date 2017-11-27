function limits=get_metric_limits(metric, varargin)

% Parse inputs
defaults = struct('symmetric', false, 'clip_vals', true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

if params.clip_vals
    pos_vals=find(metric>0);
    neg_vals=find(metric<0);
    if length(neg_vals) && length(pos_vals)
        pos_percentiles=tiedrank(metric(pos_vals))/length(pos_vals);
        neg_percentiles=tiedrank(metric(neg_vals))/length(neg_vals);
        min_clipped_val=min(metric(neg_vals(find(neg_percentiles>=.02))));
        max_clipped_val=max(metric(pos_vals(find(pos_percentiles<=.98))));
    else
        percentiles=tiedrank(metric)/length(metric);
        min_clipped_val=min(metric(find(percentiles>=.02)));
        max_clipped_val=max(metric(find(percentiles<=.98)));
    end    
    
    if params.symmetric
        limits=[-max(abs([min_clipped_val max_clipped_val])) max(abs([min_clipped_val max_clipped_val]))];
    else
        limits=[min_clipped_val max_clipped_val];
    end
else
    if params.symmetric
        limits=[-max(abs(metric)) max(abs(metric))];
    else
        limits=[min(metric) max(metric)];
    end
end