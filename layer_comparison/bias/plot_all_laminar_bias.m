function plot_all_laminar_bias(subjects, contrasts, sub_dir, varargin)

% Parse inputs
defaults = struct('data_dir','d:/pred_coding',...
    'surf_dir', 'D:/pred_coding/surf','inv_type','EBB',...
    'patch_size',0.4,'thresh_percentile',80,'roi_type','mean',...
    'whole_brain', false, 'plot_ext','',...
    'recompute_roi',false, 'ylim', []);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');
addpath('D:\pred_coding\src\matlab\analysis\layer_comparison');

plot_dir=fullfile('C:\Users\jbonai\Dropbox\meg\pred_coding\plots\layer_comparison',sub_dir);
mkdir(plot_dir);

contrast_order=[3 4 6 5 1 2];
colors=[231 138 195;102 194 165;102 194 165;231 138 195;231 138 195;102 194 165]./255;
styles={'-','-','-','--','--','--'};

all_pial_wm_diff_data={};
all_lfn_diff={};
all_depth_diff={};

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
    c_all_lfn_diff=[];
    c_all_depth_diff=[];


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
        lfn_diff=pial_lfn-mapped_wm_lfn;

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

        pial_depth_fname=fullfile(subj_surf_dir,'ds_pial_sulcal_depth.mat');
        wm_depth_fname=fullfile(subj_surf_dir,'ds_white_sulcal_depth.mat');
        if exist(pial_depth_fname,'file')~=1 || exist(wm_depth_fname,'file')~=2
            [pial_depth,HS]=compute_sulcal_depth(pial_mesh);
            mapping=dsearchn(HS.vertices,wm.vertices);
            wm_depth=sqrt(sum((wm.vertices-HS.vertices(mapping,:)).^2,2));
            save(pial_depth_fname,'pial_depth');
            save(wm_depth_fname,'wm_depth');
        else
            load(pial_depth_fname);
            load(wm_depth_fname);
        end
        mapped_wm_depth=wm_depth(pial_white_map);
        depth_diff=pial_depth-mapped_wm_depth;

        c_all_pial_wm_diff_data(end+1:end+length(mask))=pial_wm_diff_data(mask);
        c_all_lfn_diff(end+1:end+length(mask))=lfn_diff(mask);
        c_all_depth_diff(end+1:end+length(mask))=depth_diff(mask);
    end
    all_pial_wm_diff_data{c_idx}=c_all_pial_wm_diff_data;
    all_lfn_diff{c_idx}=c_all_lfn_diff;
    all_depth_diff{c_idx}=c_all_depth_diff;
end

fig=figure();
hold all;
contrast_names={};
for c_idx=1:length(contrast_order)
    contrast=contrasts(contrast_order(c_idx));
    contrast_names{c_idx}=contrast.comparison_name;
    
    depth_diff=all_depth_diff{c_idx};
    pial_wm_diff_data=all_pial_wm_diff_data{c_idx};
    
    pPoly = polyfit(depth_diff, pial_wm_diff_data,1); % Linear fit of xdata vs ydata
    linePointsX = [min(depth_diff) max(depth_diff)]; % find left and right x values
    linePointsY = polyval(pPoly,linePointsX); 
    plot (linePointsX,linePointsY,styles{c_idx},'Color',colors(c_idx,:),'LineWidth',2);
end
legend(contrast_names);
xlabel('Pial depth - WM depth','FontSize',14,'FontName','Arial');
ylabel('CV(|Pial Power D| - |WM Power D|)','FontSize',14,'FontName','Arial');
set(gca,'FontName','Arial');
saveas(fig, fullfile(plot_dir, 'depth_diff.eps'), 'epsc');


fig=figure();
hold all;
contrast_names={};
for c_idx=1:length(contrast_order)
    contrast=contrasts(contrast_order(c_idx));
    contrast_names{c_idx}=contrast.comparison_name;
    
    lfn_diff=all_lfn_diff{c_idx};
    pial_wm_diff_data=all_pial_wm_diff_data{c_idx};
    
    pPoly = polyfit(lfn_diff, pial_wm_diff_data,1); % Linear fit of xdata vs ydata
    linePointsX = [min(lfn_diff) max(lfn_diff)]; % find left and right x values
    linePointsY = polyval(pPoly,linePointsX); 
    plot (linePointsX,linePointsY,styles{c_idx},'Color',colors(c_idx,:),'LineWidth',2);
end
legend(contrast_names);
xlabel('Pial LFN - WM LFN','FontSize',14,'FontName','Arial');
ylabel('CV(|Pial Power D| - |WM Power D|)','FontSize',14,'FontName','Arial');
set(gca,'FontName','Arial');
saveas(fig, fullfile(plot_dir, 'lfn_diff.eps'), 'epsc');

rmpath('D:\pred_coding\src\matlab\analysis\layer_comparison');