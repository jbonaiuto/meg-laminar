function limits=get_pial_wm_limits(pial_metric, wm_metric, varargin)

% Parse inputs
defaults = struct('symmetric', false, 'clip_vals', true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

if params.clip_vals
    pial_pos_vals=find(pial_metric>0);
    pial_neg_vals=find(pial_metric<0);
    if length(pial_neg_vals) && length(pial_pos_vals)
        pial_pos_percentiles=tiedrank(pial_metric(pial_pos_vals))/length(pial_pos_vals);
        pial_neg_percentiles=tiedrank(pial_metric(pial_neg_vals))/length(pial_neg_vals);
        pial_min_clipped_val=min(pial_metric(pial_neg_vals(find(pial_neg_percentiles>=.02))));
        pial_max_clipped_val=max(pial_metric(pial_pos_vals(find(pial_pos_percentiles<=.98))));
    else
        pial_percentiles=tiedrank(pial_metric)/length(pial_metric);
        pial_min_clipped_val=min(pial_metric(find(pial_percentiles>=.02)));
        pial_max_clipped_val=max(pial_metric(find(pial_percentiles<=.98)));
    end

    wm_pos_vals=find(wm_metric>0);
    wm_neg_vals=find(wm_metric<0);
    if length(wm_neg_vals) && length(wm_pos_vals)
        wm_pos_percentiles=tiedrank(wm_metric(wm_pos_vals))/length(wm_pos_vals);
        wm_neg_percentiles=tiedrank(wm_metric(wm_neg_vals))/length(wm_neg_vals);
        wm_min_clipped_val=min(wm_metric(wm_neg_vals(find(wm_neg_percentiles>=.02))));
        wm_max_clipped_val=max(wm_metric(wm_pos_vals(find(wm_pos_percentiles<=.98))));
    else
        wm_percentiles=tiedrank(wm_metric)/length(wm_metric);
        wm_min_clipped_val=min(wm_metric(find(wm_percentiles>=.02)));
        wm_max_clipped_val=max(wm_metric(find(wm_percentiles<=.98)));
    end
    
    if params.symmetric
        limits=[-max(abs([pial_min_clipped_val wm_min_clipped_val pial_max_clipped_val wm_max_clipped_val])) max(abs([pial_min_clipped_val wm_min_clipped_val pial_max_clipped_val wm_max_clipped_val]))];
    else
        limits=[min([pial_min_clipped_val wm_min_clipped_val]) max([pial_max_clipped_val wm_max_clipped_val])];
    end
else
    if params.symmetric
        limits=[-max([max(pial_metric) max(wm_metric)]) max([max(pial_metric) max(wm_metric)])];
    else
        limits=[min([min(pial_metric) min(wm_metric)]) max([max(pial_metric) max(wm_metric)])];
    end
end