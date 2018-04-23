function plot_thickness_bias(subjects, contrasts, sub_dir, varargin)

% Parse inputs
defaults = struct('data_dir','d:/pred_coding',...
    'surf_dir', 'D:/pred_coding/surf','inv_type','EBB',...
    'patch_size',0.4,'thresh_percentile',80,'roi_type','mean',...
    'whole_brain', false, 'plot_ext','',...
    'recompute_roi',false, 'clim', []);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');
addpath('D:\pred_coding\src\matlab\analysis\layer_comparison');

contrast_order=[3 4 6 5 1 2];

all_pial_wm_diff_data={};
all_pial_thickness={};
    
plot_dir=fullfile('C:\Users\jbonai\Dropbox\meg\pred_coding\plots\layer_comparison');
    
    
for c_idx=1:length(contrast_order)
    contrast=contrasts(contrast_order(c_idx));
    
    thresh_type='lower';
    switch contrast.comparison_name
        case 'dots_beta_erd'
            thresh_type='upper';
        case 'dots_alpha'
            thresh_type='upper';
    end

    c_all_pial_wm_diff_data=[];
    c_all_pial_thickness=[];
    
    for s_idx=1:length(subjects)
        subj_info=subjects(s_idx);
        subj_surf_dir=fullfile(params.surf_dir,sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf');

        orig_white_mesh=fullfile(subj_surf_dir,'white.hires.deformed.surf.gii');
        white_mesh=fullfile(subj_surf_dir,'ds_white.hires.deformed.surf.gii');
        white_inflated=fullfile(subj_surf_dir,'ds_white.hires.deformed_inflated.surf.gii');
        orig_pial_mesh=fullfile(subj_surf_dir,'pial.hires.deformed.surf.gii');
        pial_mesh=fullfile(subj_surf_dir,'ds_pial.hires.deformed.surf.gii');
        pial_inflated=fullfile(subj_surf_dir,'ds_pial.hires.deformed_inflated.surf.gii');
        pial_white_map=map_pial_to_white(white_mesh, pial_mesh, 'mapType', 'link',...
            'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh);
        white_pial_map=map_white_to_pial(white_mesh, pial_mesh, 'mapType', 'link',...
            'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh);
        pial_hemisphere_map=get_hemisphere_map(pial_mesh, orig_pial_mesh);
        white_hemisphere_map=get_hemisphere_map(white_mesh, orig_white_mesh);
        pial=gifti(pial_mesh);
        wm=gifti(white_mesh);
        nvertices=size(pial.vertices,1);

        first_session=1;
        if strcmp(subj_info.subj_id,'nc')
            first_session=3;
        end
        foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
            num2str(first_session), 'grey_coreg', params.inv_type,....
            ['p' num2str(params.patch_size)], contrast.zero_event,...
            ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);

        lfn_filename=fullfile(foi_dir, sprintf('br%s_%d.mat',subj_info.subj_id, first_session));    
        grey_lfn_fname=fullfile(subj_surf_dir, 'ds_grey_lfn.mat');

        if exist(grey_lfn_fname,'file')~=2
            D=spm_eeg_load(lfn_filename);
            lfn=sqrt(sum(D.inv{1}.inverse.L.^2,1))';        
            save(grey_lfn_fname,'lfn');
        else
            load(grey_lfn_fname);
        end
        pial_lfn=lfn(nvertices+1:end);
        pial_lfn=spm_mesh_smooth(pial, pial_lfn, 20);
        wm_lfn=lfn(1:nvertices);
        wm_lfn=spm_mesh_smooth(wm, wm_lfn, 20);
        mapped_wm_lfn=wm_lfn(pial_white_map);

        pial_metric_mask=find(pial_lfn>22);
        mapped_wm_metric_mask=find(mapped_wm_lfn>19);
        metric_mask=union(pial_metric_mask,mapped_wm_metric_mask);

        hemisphere=contrast.hemisphere;
        region=contrast.region;
        if params.whole_brain
            hemisphere='';
            region='';
        end
        foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
                'grey_coreg', params.inv_type,....
                ['p' num2str(params.patch_size)], contrast.zero_event,...
                ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
        [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, contrast.comparison_name, ...
                thresh_type, pial_mesh, white_mesh, pial_inflated, white_inflated, ...
                pial_white_map, white_pial_map, lfn_filename, 'thresh_percentile',params.thresh_percentile,...
                'type',params.roi_type, 'region', region, 'hemisphere', hemisphere,...
                'pial_hemisphere_map', pial_hemisphere_map,...
                'white_hemisphere_map', white_hemisphere_map, 'recompute', params.recompute_roi);
        mask=intersect(metric_mask,mask);

        pial_wm_diff=gifti(fullfile(foi_dir,['pial-white.' contrast.comparison_name '.diff.gii']));
        pial_wm_diff_data=pial_wm_diff.cdata(:,:);
        pial_wm_diff_data=mean(pial_wm_diff_data,2)./std(pial_wm_diff_data,[],2);

        pial_thickness_fname=fullfile(subj_surf_dir,'ds_pial_thickness.mat');
        wm_thickness_fname=fullfile(subj_surf_dir,'ds_white_thickness.mat');
        if exist(pial_thickness_fname,'file')~=1 || exist(wm_thickness_fname,'file')~=2
            [pial_thickness, wm_thickness]=compute_thickness(pial_mesh, white_mesh,...
                orig_pial_mesh, orig_white_mesh);
            save(pial_thickness_fname,'pial_thickness');
            save(wm_thickness_fname,'wm_thickness');
        else
            load(pial_thickness_fname);
            load(wm_thickness_fname);
        end

        c_all_pial_wm_diff_data(end+1:end+length(mask))=pial_wm_diff_data(mask);
        c_all_pial_thickness(end+1:end+length(mask))=pial_thickness(mask);
        %all_subj(end+1:end+length(mask))=s_idx.*ones(length(mask),1);
    end
    all_pial_wm_diff_data{c_idx}=c_all_pial_wm_diff_data;
    all_pial_thickness{c_idx}=c_all_pial_thickness;
end

fig=figure();

% Red to blue colormap
r = [255 0 0]/255;       %# end
w = [255 255 255]/255;       %# middle
b = [0 0 255]/255;      %# start
c1 = zeros(128,3);
c2 = zeros(128,3);
for i=1:3
    c1(:,i) = linspace(r(i), w(i), 128);
    c2(:,i) = linspace(w(i), b(i), 128);
end
redBlue256 = [c1;c2];
colormap(redBlue256);

hold on;

xlabels={};
for c_idx=1:length(contrast_order)
    contrast=contrasts(contrast_order(c_idx));
    c_thickness=all_pial_thickness{c_idx};
    xlabels{c_idx}=contrast.comparison_name;
    
    left=c_idx-.25;
    right=c_idx+.25;
    data_x=left+(right-left).*rand([length(c_thickness) 1]);
    scatter(data_x,c_thickness,25,all_pial_wm_diff_data{c_idx},'.');       
end
if length(params.clim)>0
    set(gca,'clim',params.clim);
end
set(gca,'xtick',[1:length(contrast_order)]);
set(gca,'xticklabel', xlabels);
ylabel('Cortical thickness (mm)', 'FontName','arial','FontSize',14);
hcb=colorbar();
set(get(hcb,'Title'),'String','CV(|Pial Power D| - |WM Power D|)');
set(gca,'FontName','Arial');
saveas(fig, fullfile(plot_dir, 'thickness_bias.eps'), 'epsc');
saveas(fig, fullfile(plot_dir, 'thickness_bias.png'), 'png');   

fig=figure();
for c_idx=1:length(contrast_order)
    subplot(length(contrast_order),1,c_idx);
    hold all
    contrast=contrasts(contrast_order(c_idx));
    c_thickness=all_pial_thickness{c_idx};
    
    dx = 0.1*range(c_thickness) ;
    xLim = [min(c_thickness)-dx, max(c_thickness)+dx];
    xrange = xLim(1):0.01*dx: xLim(2);
    
    classified_pial=find(all_pial_wm_diff_data{c_idx}>0);
    pial_px = ksdensity(c_thickness(classified_pial),xrange)*(length(classified_pial)/length(all_pial_wm_diff_data{c_idx}));
    plot(xrange,pial_px,'Color',[0 0 255]/255,'LineWidth',2);
    
    classified_wm=find(all_pial_wm_diff_data{c_idx}<0);
    white_px = ksdensity(c_thickness(classified_wm),xrange)*(length(classified_wm)/length(all_pial_wm_diff_data{c_idx}));
    plot(xrange,white_px,'Color',[255 0 0]/255,'LineWidth',2);
    
    xlim([1.5 6.5]);
    ylabel(contrast.comparison_name,'FontName', 'arial', 'FontSize', 14);
end
xlabel('Cortical thickness (mm)', 'FontName', 'arial', 'FontSize', 14);
saveas(fig, fullfile(plot_dir, 'thickness_bias_hist.eps'), 'epsc');
saveas(fig, fullfile(plot_dir, 'thickness_bias_hist.png'), 'png');   


rmpath('D:\pred_coding\src\matlab\analysis\layer_comparison');