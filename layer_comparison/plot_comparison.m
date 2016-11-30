function [pial_data,white_data,pial_white_data,pial_white_diff_data]=plot_comparison(subj_info, foi_dir, comparison_name, view, varargin)

% Parse inputs
defaults = struct('subjects_dir',fullfile('/usr','local','freesurfer','subjects'), 'output_file','','output_format','png','threshold',0,'hist_ylim',0,'hist_xlim',[-25 25], 'white_pial_map', [], 'plot', true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

pial_surf=fullfile(params.subjects_dir, [subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed_inflated.surf.gii');
wm_surf=fullfile(params.subjects_dir, [subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed_inflated.surf.gii');

pial=gifti(pial_surf);
wm=gifti(wm_surf);

pial_comparison=fullfile(foi_dir,['pial.' comparison_name '.diff.shape.gii']);
wm_comparison=fullfile(foi_dir,['white.' comparison_name '.diff.shape.gii']);
pial_wm_comparison=fullfile(foi_dir,['pial-white.' comparison_name '.t.shape.gii']);
pial_wm_diff=fullfile(foi_dir,['pial-white.' comparison_name '.diff.shape.gii']);

pial_metric=gifti(pial_comparison);
pial_pos_vals=find(pial_metric.cdata(:)>0);
pial_neg_vals=find(pial_metric.cdata(:)<0);
if length(pial_neg_vals) && length(pial_pos_vals)
    pial_pos_percentiles=tiedrank(pial_metric.cdata(pial_pos_vals))/length(pial_pos_vals);
    pial_neg_percentiles=tiedrank(pial_metric.cdata(pial_neg_vals))/length(pial_neg_vals);
    pial_min_clipped_val=min(pial_metric.cdata(pial_neg_vals(find(pial_neg_percentiles>=.02))));
    pial_max_clipped_val=max(pial_metric.cdata(pial_pos_vals(find(pial_pos_percentiles<=.98))));
else
    pial_percentiles=tiedrank(pial_metric.cdata(:))/length(pial_metric.cdata(:));
    pial_min_clipped_val=min(pial_metric.cdata(find(pial_percentiles>=.02)));
    pial_max_clipped_val=max(pial_metric.cdata(find(pial_percentiles<=.98)));
end

wm_metric=gifti(wm_comparison);
wm_pos_vals=find(wm_metric.cdata(:)>0);
wm_neg_vals=find(wm_metric.cdata(:)<0);
if length(wm_neg_vals) && length(wm_pos_vals)
    wm_pos_percentiles=tiedrank(wm_metric.cdata(wm_pos_vals))/length(wm_pos_vals);
    wm_neg_percentiles=tiedrank(wm_metric.cdata(wm_neg_vals))/length(wm_neg_vals);
    wm_min_clipped_val=min(wm_metric.cdata(wm_neg_vals(find(wm_neg_percentiles>=.02))));
    wm_max_clipped_val=max(wm_metric.cdata(wm_pos_vals(find(wm_pos_percentiles<=.98))));
else
    wm_percentiles=tiedrank(wm_metric.cdata(:))/length(wm_metric.cdata(:));
    wm_min_clipped_val=min(wm_metric.cdata(find(wm_percentiles>=.02)));
    wm_max_clipped_val=max(wm_metric.cdata(find(wm_percentiles<=.98)));
end

limits=[min([pial_min_clipped_val wm_min_clipped_val]) max([pial_max_clipped_val wm_max_clipped_val])];
%limits=[min([pial_min_clipped_val wm_min_clipped_val]) 1];
ax=0;
if params.plot || length(params.output_file)>0
    fig=figure('Position',[1 1 1800 400],'PaperUnits','points','PaperPosition',[1 1 900 200],'PaperPositionMode','manual');
    ax=subplot(1,4,1);
end
pial_data=plot_surface_metric(subj_info, pial, pial_metric, view, 'ax', ax, 'limits', limits, 'threshold', params.threshold, 'plot', params.plot);

if params.plot || length(params.output_file)>0
    freezeColors;
    ax=subplot(1,4,2);
end
white_data=plot_surface_metric(subj_info, wm, wm_metric, view, 'ax', ax, 'limits', limits, 'threshold', params.threshold, 'plot', params.plot);

if params.plot || length(params.output_file)>0
    freezeColors;
    ax=subplot(1,4,3);
end

mask=[];
if abs(params.threshold)>0
    if length(params.white_pial_map)==0
        params.white_pial_map=map_white_to_pial(wm, pial);
    end
    mapped_white_metric=wm_metric.cdata(params.white_pial_map);
    if params.threshold>0
        mask=union(find(pial_metric.cdata(:)>params.threshold),find(mapped_white_metric>params.threshold));
    else
        mask=union(find(pial_metric.cdata(:)<params.threshold),find(mapped_white_metric<params.threshold));
    end
end
pial_white_data=plot_surface_metric(subj_info, pial, gifti(pial_wm_comparison), view, 'ax', ax, 'mask', mask, 'plot', params.plot);
pial_white_diff_surf=gifti(pial_wm_diff);
pial_white_diff_data=pial_white_diff_surf.cdata(:);
if length(mask)>0
    pial_white_diff_data=pial_white_diff_data(mask);
end

    
if params.plot || length(params.output_file)>0
    freezeColors;
    ax=subplot(1,4,4);
    %[f,x]=hist(pial_white_data,length(pial_white_data)/60);
    bin_width=0.5;
    [f,x]=hist(pial_white_data,[min(pial_white_data):bin_width:max(pial_white_data)]);      
    bar(x,f/sum(f).*100.0,'b');
    xlim(params.hist_xlim);
    if params.hist_ylim>0
        ylim([0 params.hist_ylim]);
    end
    xlabel('Pial-White');
    ylabel('Percent vertices');

    if length(params.output_file)>0
        saveas(fig, params.output_file, params.output_format);
    end
end
