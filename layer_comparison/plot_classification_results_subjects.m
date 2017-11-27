function plot_classification_results_subjects(subjects, contrast, varargin)

% Parse inputs
defaults = struct('data_dir','d:/pred_coding',...
    'surf_dir', 'D:/pred_coding/surf','inv_type','EBB',...
    'patch_size',0.4,'recompute_roi',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

subject_whole_brain_tvals=zeros(1,length(subjects));
subject_func_roi_tvals=zeros(1,length(subjects));
subject_anatfunc_roi_tvals=zeros(1,length(subjects));

thresh_type='lower';
switch contrast.comparison_name
    case 'dots_beta_erd'
        thresh_type='upper';
    case 'dots_alpha'
        thresh_type='upper';
end
subj_dofs=zeros(1,length(subjects));

for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
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
    
    good_sessions=[1:length(subj_info.sessions)];
    good_sessions=setdiff(good_sessions,subj_info.exclude_sessions(contrast.comparison_name));

    foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
            num2str(good_sessions(1)), 'grey_coreg', params.inv_type,....
            ['p' num2str(params.patch_size)], contrast.zero_event,...
            ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
    lfn_filename=fullfile(foi_dir, sprintf('br%s_%d.mat',subj_info.subj_id, good_sessions(1)));    
    
    foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
            'grey_coreg', params.inv_type,....
            ['p' num2str(params.patch_size)], contrast.zero_event,...
            ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
    % Get mean pial-wm in ROI
    pial_wm_diff=gifti(fullfile(foi_dir,['pial-white.' contrast.comparison_name '.diff.gii']));
            
    subj_dofs(subj_idx)=size(pial_wm_diff.cdata(:,:),2)-1;
    
    % whole brain
    [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, contrast.comparison_name, ...
        thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
        pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',0,...
        'type','mean', 'region', '', 'hemisphere', '',...
        'pial_hemisphere_map', pial_hemisphere_map,...
        'white_hemisphere_map', white_hemisphere_map, 'recompute', params.recompute_roi);            
    pial_wm_roi_diff=mean(pial_wm_diff.cdata(mask,:));
    [tstat,p]=ttest_corrected(pial_wm_roi_diff','correction',25*var(pial_wm_roi_diff));
    disp(sprintf('%s, whole=%.3f',subj_info.subj_id,tstat));
    subject_whole_brain_tvals(subj_idx)=tstat;
    
    % func ROI
    [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, contrast.comparison_name, ...
        thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
        pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',80,...
        'type','mean', 'region', '', 'hemisphere', '',...
        'pial_hemisphere_map', pial_hemisphere_map,...
        'white_hemisphere_map', white_hemisphere_map, 'recompute', params.recompute_roi);            
    pial_wm_roi_diff=mean(pial_wm_diff.cdata(mask,:));
    [tstat,p]=ttest_corrected(pial_wm_roi_diff','correction',25*var(pial_wm_roi_diff));
    disp(sprintf('%s, func=%.3f',subj_info.subj_id,tstat));
    subject_func_roi_tvals(subj_idx)=tstat;
    
    % anat-func ROI
    [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, contrast.comparison_name, ...
        thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
        pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',80,...
        'type','mean', 'region', contrast.region, 'hemisphere', contrast.hemisphere,...
        'pial_hemisphere_map', pial_hemisphere_map,...
        'white_hemisphere_map', white_hemisphere_map, 'recompute', params.recompute_roi);            
    pial_wm_roi_diff=mean(pial_wm_diff.cdata(mask,:));
    [tstat,p]=ttest_corrected(pial_wm_roi_diff','correction',25*var(pial_wm_roi_diff));
    disp(sprintf('%s, anat-func=%.3f',subj_info.subj_id,tstat));
    subject_anatfunc_roi_tvals(subj_idx)=tstat;    
end

fig=figure('Position',[1 1 1200 600],'PaperUnits','points',...
    'PaperPosition',[1 1 600 300],'PaperPositionMode','manual');
hold on;
bar_width=0.1;
gap_width=0.05;
subj_width=3*bar_width+(3-1)*gap_width;

subj_ids={};
alpha=1.0-(0.05/2);
t_thresh=tinv(alpha, mean(subj_dofs));
plot([1-.5 length(subjects)+.5],[t_thresh t_thresh],'k--','LineWidth',2);
plot([1-.5 length(subjects)+.5],[-t_thresh -t_thresh],'k--','LineWidth',2);

    
for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    subj_ids{subj_idx}=num2str(subj_idx);
    left=subj_idx-.5*subj_width;
        
    center=left+.5*bar_width+(1-1)*(bar_width+gap_width);
    tval=subject_whole_brain_tvals(subj_idx);
    face_color=[0 39 102]./255.0;
    if tval<0
        face_color=[130 0 0]./255.0;
    end
    bar(center, tval, bar_width, 'FaceColor', face_color,'EdgeColor','none');
    
    center=left+.5*bar_width+(2-1)*(bar_width+gap_width);
    tval=subject_func_roi_tvals(subj_idx);
    face_color=[0 64 168]./255.0;
    if tval<0
        face_color=[200 0 0]./255.0;
    end
    bar(center, tval, bar_width, 'FaceColor', face_color,'EdgeColor','none');
    
    center=left+.5*bar_width+(3-1)*(bar_width+gap_width);
    tval=subject_anatfunc_roi_tvals(subj_idx);
    face_color=[0 97 255]./255.0;
    if tval<0
        face_color=[255 0 0]./255.0;
    end
    bar(center, tval, bar_width, 'FaceColor', face_color,'EdgeColor','none');            
end
set(gca,'XTick',1:length(subjects));
set(gca,'XTickLabel',subj_ids);
xlabel('Participant','Fontsize',24,'Fontname','Arial');
ylabel('Pial-White ROI t-statistic','Fontsize',24,'Fontname','Arial');
xlim([.5 length(subjects)+.5]);
yl=ylim;
ylim([-max(abs(yl)) max(abs(yl))]);
xt=get(gca,'XTick');
set(gca,'FontSize',20);
set(gca,'Fontname','Arial');

plot_dir=fullfile('C:\Users\jbonai\Dropbox\meg\pred_coding\plots\layer_comparison',...
        contrast.comparison_name);
mkdir(plot_dir);
saveas(fig, fullfile(plot_dir, sprintf('%s_subjects.png', contrast.comparison_name)), 'png');
saveas(fig, fullfile(plot_dir, sprintf('%s_subjects.eps', contrast.comparison_name)), 'eps');
saveas(fig, fullfile(plot_dir, sprintf('%s_subjects.fig', contrast.comparison_name)), 'fig');
