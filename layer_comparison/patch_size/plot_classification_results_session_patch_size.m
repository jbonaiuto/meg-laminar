function plot_classification_results_session_patch_size(subj_info, session_num, contrasts, patch_sizes, varargin)

% Parse inputs
defaults = struct('data_dir','d:/pred_coding',...
    'surf_dir', 'D:/pred_coding/surf','inv_type','EBB',...
    'thresh_percentile',80,'roi_type','mean',...
    'whole_brain', false, 'plot_ext','',...
    'recompute_roi',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

addpath('D:\pred_coding\src\matlab\analysis\layer_comparison');

spm('defaults','eeg');

orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');
white_inflated=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed_inflated.surf.gii');
orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');
pial_inflated=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed_inflated.surf.gii');
pial_white_map=map_pial_to_white(white_mesh, pial_mesh, 'mapType', 'link',...
    'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh);
white_pial_map=map_white_to_pial(white_mesh, pial_mesh, 'mapType', 'link',...
    'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh);
pial_hemisphere_map=get_hemisphere_map(pial_mesh, orig_pial_mesh);
white_hemisphere_map=get_hemisphere_map(white_mesh, orig_white_mesh);
    
contrast_patch_sizes=zeros(length(contrasts), length(patch_sizes));

for contrast_idx=1:length(contrasts)
    contrast=contrasts(contrast_idx);

    view='back';
    thresh_type='lower';
    switch contrast.comparison_name
        case 'dots_beta_erd'
            view='topdown';
            thresh_type='upper';
        case 'resp_mrgs'
            view='topdown';
        case 'resp_beta_rebound'
            view='topdown';
        case 'dots_alpha'
            thresh_type='upper';
    end
    
    for patch_size_idx=1:length(patch_sizes)
        patch_size=patch_sizes(patch_size_idx);        
    
        foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
                num2str(session_num), 'grey_coreg', params.inv_type,....
                ['p' num2str(patch_size)], contrast.zero_event,...
                ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
        lfn_filename=fullfile(foi_dir, sprintf('br%s_%d.mat',subj_info.subj_id, 1));
    
        plot_dir=fullfile('C:\Users\jbonai\Dropbox\meg\pred_coding\plots\layer_comparison',...
            subj_info.subj_id, num2str(session_num), contrast.comparison_name, ['p' num2str(patch_size)]);
        mkdir(plot_dir);
        
        unthresh_out_file=fullfile(plot_dir,sprintf('%s_%sunthresholded.png',contrast.comparison_name,params.plot_ext));
        [diff_limits, pial_wm_diff_limits, t_limits, pial_wm_t_limits]=plot_comparison(subj_info,...
            lfn_filename,foi_dir,contrast.comparison_name,view,...
            'surf_dir',params.surf_dir,'threshold',[]);
        unthresh_fig=gcf;
        unthresh_children=get(unthresh_fig,'Children');
        unthresh_t_ax=unthresh_children(8);
        unthresh_hist_t_ylim=get(unthresh_t_ax,'ylim');
    
        hemisphere=contrast.hemisphere;
        region=contrast.region;
        if params.whole_brain
            hemisphere='';
            region='';
        end
        [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, contrast.comparison_name, ...
            thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
            pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',params.thresh_percentile,...
            'type',params.roi_type, 'region', region, 'hemisphere', hemisphere,...
            'pial_hemisphere_map', pial_hemisphere_map,...
            'white_hemisphere_map', white_hemisphere_map,'recompute',params.recompute_roi);
        
        thresh_out_file=fullfile(plot_dir,sprintf('%s_%sthresholded.png',contrast.comparison_name,params.plot_ext));
        plot_comparison(subj_info,lfn_filename,foi_dir,...
            contrast.comparison_name,view,'surf_dir',params.surf_dir,...
            'pial_mask',pial_mask,'wm_mask',wm_mask,'mask',mask,...
            'diff_limits', diff_limits, 't_limits', t_limits,...
            'pial_wm_diff_limits', pial_wm_diff_limits,...
            'pial_wm_t_limits', pial_wm_t_limits);
        thresh_fig=gcf;
        thresh_children=get(thresh_fig,'Children');
        thresh_t_ax=thresh_children(8);
        thresh_hist_t_ylim=get(thresh_t_ax,'ylim');
        hist_t_ylim=[min([unthresh_hist_t_ylim thresh_hist_t_ylim]) max([unthresh_hist_t_ylim thresh_hist_t_ylim])];
        set(unthresh_t_ax,'ylim',hist_t_ylim);
        set(thresh_t_ax,'ylim',hist_t_ylim);

        saveas(unthresh_fig, unthresh_out_file, 'png');
        saveas(thresh_fig, thresh_out_file, 'png');

        % Get mean pial-wm in ROI
        pial_wm_diff=gifti(fullfile(foi_dir,['pial-white.' contrast.comparison_name '.diff.gii']));
        pial_wm_roi_diff=mean(pial_wm_diff.cdata(mask,:));
        % Perform ROI t-stat
        [H,pvals,ci,STATS]=ttest(pial_wm_roi_diff');
        tstat=STATS.tstat;
        disp(sprintf('%s: %.2f=%.3f',contrast.comparison_name,patch_size,tstat));
        contrast_patch_sizes(contrast_idx,patch_size_idx)=tstat;        
    end
    close all;
end

fig=figure();
hold on;
bar_width=0.1;
gap_width=0.05;
contrast_names={};
for contrast_idx=1:length(contrasts)
    contrast=contrasts(contrast_idx);
    contrast_names{contrast_idx}=contrast.comparison_name;
    contrast_width=length(patch_sizes)*bar_width+(length(patch_sizes)-1)*gap_width;
    left=contrast_idx-.5*contrast_width;
    for patch_size_idx=1:length(patch_sizes)
        center=left+.5*bar_width+(patch_size_idx-1)*(bar_width+gap_width);
        tval=contrast_patch_sizes(contrast_idx,patch_size_idx);
        color='b';
        if tval<0
            color='r';
        end
        bar(center, tval, bar_width, 'FaceColor', color,'EdgeColor','none');
    end
end
set(gca,'XTick',1:length(contrasts));
set(gca,'XTickLabel',contrast_names);
yl=ylim;
ylim([-max(abs(yl)) max(abs(yl))]);
ylabel('Pial-White ROI t-statistic');
plot_dir=fullfile('C:\Users\jbonai\Dropbox\meg\pred_coding\plots\layer_comparison',...
        subj_info.subj_id, num2str(session_num), 'patch_size');
mkdir(plot_dir);
saveas(fig, fullfile(plot_dir, sprintf('patch_size_%ssessions.png', params.plot_ext)), 'png');

rmpath('D:\pred_coding\src\matlab\analysis\layer_comparison');