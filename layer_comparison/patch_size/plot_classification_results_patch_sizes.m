function plot_classification_results_patch_sizes(subj_info, contrasts, roi_type, varargin)

% Parse inputs
defaults = struct('data_dir','d:/pred_coding',...
    'surf_dir', 'D:/pred_coding/surf','inv_type','EBB',...
    'patch_size',0.4,'recompute_roi',false, 'whole_brain', false,...
    'thresh_percentile', 80);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');
addpath('D:\pred_coding\src\matlab\analysis\layer_comparison');

patch_sizes=[2.5 5 10 20];
tvals=zeros(length(contrasts),length(patch_sizes));
mean_f=zeros(length(contrasts),length(patch_sizes));
mean_pial_f=zeros(length(contrasts),length(patch_sizes));
mean_white_f=zeros(length(contrasts),length(patch_sizes));
contrast_order=[3 4 6 5 1 2];

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
    
contrast_dofs=zeros(length(contrasts),length(patch_sizes));

for c_idx=1:length(contrast_order)
    contrast=contrasts(contrast_order(c_idx));
    
    thresh_type='lower';
    switch contrast.comparison_name
        case 'dots_beta_erd'
            thresh_type='upper';
        case 'dots_alpha'
            thresh_type='upper';
    end
    region='';
    hemisphere='';
    if ~params.whole_brain
        region=contrast.region;
        hemisphere=contrast.hemisphere;
    end    
    
    for p_idx=1:length(patch_sizes)
        patch_size=patch_sizes(p_idx);
        
        foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
                num2str(subj_info.sessions(1)), 'grey_coreg', params.inv_type,....
                ['p' num2str(patch_size)], contrast.zero_event,...
                ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
        lfn_filename=fullfile(foi_dir, sprintf('br%s_%d.mat',subj_info.subj_id, subj_info.sessions(1)));    

        foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
                'grey_coreg', params.inv_type,....
                ['p' num2str(patch_size)], contrast.zero_event,...
                ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
        % Get mean pial-wm in ROI
        pial_wm_diff=gifti(fullfile(foi_dir,['pial-white.' contrast.comparison_name '.diff.gii']));

        contrast_dofs(c_idx, p_idx)=size(pial_wm_diff.cdata(:,:),2)-1;

        [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, contrast.comparison_name, ...
            thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
            pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',params.thresh_percentile,...
            'type','mean', 'region', region, 'hemisphere', hemisphere,...
            'pial_hemisphere_map', pial_hemisphere_map,...
            'white_hemisphere_map', white_hemisphere_map, 'recompute', params.recompute_roi);            
        pial_wm_roi_diff=mean(pial_wm_diff.cdata(mask,:));
        [tstat,p]=ttest_corrected(pial_wm_roi_diff','correction',25*var(pial_wm_roi_diff));
        disp(sprintf('%s, p=%.2f, t=%.3f',contrast.comparison_name,patch_size,tstat));
        tvals(c_idx,p_idx)=tstat;        
        
        f=zeros(1,length(subj_info.sessions));
        for session_idx=1:length(subj_info.sessions)
            session_num=subj_info.sessions(session_idx);
            foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
                num2str(session_num), 'grey_coreg', params.inv_type,....
                ['p' num2str(patch_size)], contrast.zero_event,...
                ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
            load(fullfile(foi_dir, sprintf('br%s_%d.mat',subj_info.subj_id, session_num)));
            f(session_idx)=D.other.inv{1}.inverse.F;
        end
        mean_f(c_idx,p_idx)=mean(f);
        
        f=zeros(1,length(subj_info.sessions));
        for session_idx=1:length(subj_info.sessions)
            session_num=subj_info.sessions(session_idx);
            foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
                num2str(session_num), 'pial_coreg', params.inv_type,....
                ['p' num2str(patch_size)], contrast.zero_event,...
                ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
            load(fullfile(foi_dir, sprintf('br%s_%d.mat',subj_info.subj_id, session_num)));
            f(session_idx)=D.other.inv{1}.inverse.F;
        end
        mean_pial_f(c_idx,p_idx)=mean(f);
        
        f=zeros(1,length(subj_info.sessions));
        for session_idx=1:length(subj_info.sessions)
            session_num=subj_info.sessions(session_idx);
            foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
                num2str(session_num), 'white_coreg', params.inv_type,....
                ['p' num2str(patch_size)], contrast.zero_event,...
                ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
            load(fullfile(foi_dir, sprintf('br%s_%d.mat',subj_info.subj_id, session_num)));
            f(session_idx)=D.other.inv{1}.inverse.F;
        end
        mean_white_f(c_idx,p_idx)=mean(f);
    end
end

plot_dir=fullfile('C:\Users\jbonai\Dropbox\meg\pred_coding\plots\layer_comparison\patch_size');
mkdir(plot_dir);


fig=figure('Position',[1 1 1200 600],'PaperUnits','points',...
    'PaperPosition',[1 1 600 300],'PaperPositionMode','manual');
hold on;
bar_width=0.1;
gap_width=0.05;
contrast_width=length(patch_sizes)*bar_width+(length(patch_sizes)-1)*gap_width;
contrast_names={};
for c_idx=1:length(contrast_order)
    contrast=contrasts(contrast_order(c_idx));
    contrast_names{c_idx}=contrast.comparison_name;
    left=c_idx-.5*contrast_width;
    centers=left+.5*bar_width+([1:length(patch_sizes)]-1)*(bar_width+gap_width);
    contrast_mean_f=mean_f(c_idx,:);
    plot(centers,contrast_mean_f-mean(contrast_mean_f),'k','LineWidth',2);
    contrast_mean_f=mean_pial_f(c_idx,:);
    plot(centers,contrast_mean_f-mean(contrast_mean_f),'b','LineWidth',2);
    contrast_mean_f=mean_white_f(c_idx,:);
    plot(centers,contrast_mean_f-mean(contrast_mean_f),'r','LineWidth',2);
end
legend('Grey','Pial','White');
set(gca,'XTick',1:length(contrasts));
set(gca,'XTickLabel',contrast_names);
ylabel('D F','Fontsize',24,'Fontname','Arial');
xlim([.5 length(contrasts)+.5]);
yl=ylim;
ylim([-max(abs(yl)) max(abs(yl))]);
set(gca,'FontSize',20);
set(gca,'Fontname','Arial');

saveas(fig, fullfile(plot_dir, sprintf('%s_%s_patch_size_F.png', subj_info.subj_id, roi_type)), 'png');
saveas(fig, fullfile(plot_dir, sprintf('%s_%s_patch_size_F.eps', subj_info.subj_id, roi_type)), 'eps');
saveas(fig, fullfile(plot_dir, sprintf('%s_%s_patch_size_F.fig', subj_info.subj_id, roi_type)), 'fig');


fig=figure('Position',[1 1 1200 600],'PaperUnits','points',...
    'PaperPosition',[1 1 600 300],'PaperPositionMode','manual');
hold on;
bar_width=0.1;
gap_width=0.05;
contrast_width=length(patch_sizes)*bar_width+(length(patch_sizes)-1)*gap_width;

contrast_names={};
alpha=1.0-(0.05/2);
t_thresh=tinv(alpha, mean(contrast_dofs(:)));
plot([1-.5 length(contrasts)+.5],[t_thresh t_thresh],'k--','LineWidth',2);
plot([1-.5 length(contrasts)+.5],[-t_thresh -t_thresh],'k--','LineWidth',2);

    
colors=[178,24,43; 214,96,77; 244,165,130; 253,219,199; 209,229,240; 146,197,222; 67,147,195; 33,102,172]./255;

for c_idx=1:length(contrast_order)
    contrast=contrasts(contrast_order(c_idx));
    contrast_names{c_idx}=contrast.comparison_name;
    left=c_idx-.5*contrast_width;
        
    for p_idx=1:length(patch_sizes)
        center=left+.5*bar_width+(p_idx-1)*(bar_width+gap_width);
        tval=tvals(c_idx,p_idx);
        face_color=colors(p_idx+4,:);
        if tval<0
            face_color=colors(p_idx,:);
        end
        bar(center, tval, bar_width, 'FaceColor', face_color,'EdgeColor','none');
    end    
end
set(gca,'XTick',1:length(contrasts));
set(gca,'XTickLabel',contrast_names);
ylabel('Pial-White ROI t-statistic','Fontsize',24,'Fontname','Arial');
xlim([.5 length(contrasts)+.5]);
yl=ylim;
ylim([-max(abs(yl)) max(abs(yl))]);
set(gca,'FontSize',20);
set(gca,'Fontname','Arial');

saveas(fig, fullfile(plot_dir, sprintf('%s_%s_patch_size.png', subj_info.subj_id, roi_type)), 'png');
saveas(fig, fullfile(plot_dir, sprintf('%s_%s_patch_size.eps', subj_info.subj_id, roi_type)), 'eps');
saveas(fig, fullfile(plot_dir, sprintf('%s_%s_patch_size.fig', subj_info.subj_id, roi_type)), 'fig');

rmpath('D:\pred_coding\src\matlab\analysis\layer_comparison');