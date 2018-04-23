function [diff_limits, pial_wm_diff_limits, t_limits, pial_wm_t_limits]=plot_comparison(subj_info, lfn_filename, foi_dir, comparison_name, view, varargin)

% Parse inputs
defaults = struct('data_dir','d:/pred_coding','surf_dir', '', 'mri_dir', '',...
    'inv_type','EBB','patch_size',0.4,'output_file','','output_format','png',...
    'threshold',0,'thresh_type','lower','pial_mask',[],'wm_mask',[],'mask',[],...
    't_limits', [], 'diff_limits', [], 'pial_wm_t_limits', [], 'pial_wm_diff_limits', [],...
    'hist_ylim',0,'hist_xlim',[-25 25]);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end
if length(params.surf_dir)==0
    params.surf_dir=fullfile(params.data_dir,'surf');
end
if length(params.mri_dir)==0
    params.mri_dir=fullfile(params.data_dir,'mri');
end

% Load surface files
orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');
white_inflated_mesh=gifti(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed_inflated.surf.gii'));
orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');
pial_inflated_mesh=gifti(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed_inflated.surf.gii'));
    
pial_white_map=map_pial_to_white(white_mesh, pial_mesh, 'mapType', 'link',...
    'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh);

pial_t=gifti(fullfile(foi_dir,['pial.' comparison_name '.t.gii']));
pial_t_data=pial_t.cdata(:);
wm_t=gifti(fullfile(foi_dir,['white.' comparison_name '.t.gii']));
wm_t_data=wm_t.cdata(:);
pial_wm_t=gifti(fullfile(foi_dir,['pial-white.' comparison_name '.t.gii']));
pial_wm_t_data=pial_wm_t.cdata(:);

pial_diff=gifti(fullfile(foi_dir,['pial.' comparison_name '.diff.gii']));
pial_diff_data=pial_diff.cdata(:,:);
wm_diff=gifti(fullfile(foi_dir,['white.' comparison_name '.diff.gii']));
wm_diff_data=wm_diff.cdata(:,:);
pial_wm_diff=gifti(fullfile(foi_dir,['pial-white.' comparison_name '.diff.gii']));
pial_wm_diff_data=pial_wm_diff.cdata(:,:);

if length(params.pial_mask)==0 || length(params.wm_mask)==0 || length(params.mask)==0
    [params.pial_mask,params.wm_mask,params.mask]=get_pial_wm_mask(pial_t_data, wm_t_data,...
        params.threshold, pial_white_map, 'thresh_type', params.thresh_type);
end

pial=gifti(pial_mesh);
wm=gifti(white_mesh);
nvertices=size(pial.vertices,1);

D=spm_eeg_load(lfn_filename);
pial_lfn=sqrt(sum(D.inv{1}.inverse.L(:,nvertices+1:end).^2,1))';
pial_lfn=spm_mesh_smooth(pial, pial_lfn, 20);

wm_lfn=sqrt(sum(D.inv{1}.inverse.L(:,1:nvertices).^2,1))';
wm_lfn=spm_mesh_smooth(wm, wm_lfn, 20);

mapped_wm_lfn=wm_lfn(pial_white_map);
pial_metric_mask=find(pial_lfn>22);
wm_metric_mask=find(wm_lfn>19);
mapped_wm_metric_mask=find(mapped_wm_lfn>19);

if length(params.mask)
    params.pial_mask=intersect(params.pial_mask, pial_metric_mask);
    params.wm_mask=intersect(params.wm_mask, wm_metric_mask);
    params.mask=intersect(params.mask,union(pial_metric_mask,mapped_wm_metric_mask));
else
    params.pial_mask=pial_metric_mask;
    params.wm_mask=wm_metric_mask;
    params.mask=union(pial_metric_mask,mapped_wm_metric_mask);
end

if length(params.t_limits)==0
    params.t_limits=get_pial_wm_limits(pial_t_data(params.pial_mask), ...
        wm_t_data(params.wm_mask), 'clip_vals', true, 'symmetric', true);
end
t_limits=params.t_limits;
if length(params.pial_wm_t_limits)==0
    params.pial_wm_t_limits=get_metric_limits(pial_wm_t_data(params.mask),...
        'clip_vals', false, 'symmetric', true);
end
pial_wm_t_limits=params.pial_wm_t_limits;
if length(params.diff_limits)==0
    params.diff_limits=get_pial_wm_limits(mean(pial_diff_data(params.pial_mask,:),2), ...
        mean(wm_diff_data(params.wm_mask,:),2), 'clip_vals', true, 'symmetric', true);                    
end
if length(params.pial_wm_diff_limits)==0
    params.pial_wm_diff_limits=get_metric_limits(mean(pial_wm_diff_data(params.mask,:),2),...
        'clip_vals', true, 'symmetric', true);
end
diff_limits=params.diff_limits;
pial_wm_diff_limits=params.pial_wm_diff_limits;

%fig=figure('Position',[1 1 1900 600],'PaperUnits','points',...
%    'PaperPosition',[1 1 900 300],'PaperPositionMode','manual');
fig=figure('Position',[1 1 1900 400],'PaperUnits','points',...
    'PaperPosition',[1 1 900 200],'PaperPositionMode','manual');
% ax=subplot(2,4,1);
% %ax=0;
% [ax,~]=plot_surface_metric(pial_inflated_mesh, pial_t_data, 'ax', ax,...
%     'limits', t_limits, 'mask', params.pial_mask);
% set(ax,'CameraViewAngle',6.028);
% set(ax,'CameraUpVector',subj_info.camera_up_vector(view));
% set(ax,'CameraPosition',subj_info.camera_position(view));
% freezeColors;
% 
% ax=subplot(2,4,2);
% %ax=0;
% [ax,~]=plot_surface_metric(white_inflated_mesh, wm_t_data, 'ax', ax,...
%     'limits', t_limits, 'mask', params.wm_mask);
% set(ax,'CameraViewAngle',6.028);
% set(ax,'CameraUpVector',subj_info.camera_up_vector(view));
% set(ax,'CameraPosition',subj_info.camera_position(view));
% freezeColors;
% 
% ax=subplot(2,4,3);
% % ax=0;
% [ax,~]=plot_surface_metric(pial_inflated_mesh, pial_wm_t_data, 'ax', ax,...
%     'limits', pial_wm_t_limits, 'mask', params.mask);
% set(ax,'CameraViewAngle',6.028);
% set(ax,'CameraUpVector',subj_info.camera_up_vector(view));
% set(ax,'CameraPosition',subj_info.camera_position(view));
% freezeColors;

%ax=subplot(2,4,4);
ax=subplot(1,4,4);
%figure();
nbins=100.0;
masked_pial_wm_t=pial_wm_t_data;
if length(params.mask)
    masked_pial_wm_t=pial_wm_t_data(params.mask);
end
bin_width=(max(masked_pial_wm_t)-min(masked_pial_wm_t))/nbins;
[f,x]=hist(masked_pial_wm_t,[min(masked_pial_wm_t):bin_width:max(masked_pial_wm_t)]);      
bar(x,f/sum(f).*100.0,'b');
xlim(pial_wm_t_limits);
%if params.hist_ylim>0
%    ylim([0 params.hist_ylim]);
%end
xlabel('Pial-White');
ylabel('Percent vertices');

%skekurtest(masked_pial_wm_t)
%[p,h,stats] = signrank(masked_pial_wm_t)

%ax=subplot(2,4,5);
ax=subplot(1,4,1);
%ax=0;
[ax,~]=plot_surface_metric(pial_inflated_mesh, mean(pial_diff_data,2), 'ax', ax,...
    'mask', params.pial_mask,'clip_vals', false, 'limits',diff_limits);
set(ax,'CameraViewAngle',6.028);
set(ax,'CameraUpVector',subj_info.camera_up_vector(view));
set(ax,'CameraPosition',subj_info.camera_position(view));
freezeColors;

%ax=subplot(2,4,6);
ax=subplot(1,4,2);
%ax=0;
[ax,~]=plot_surface_metric(white_inflated_mesh, mean(wm_diff_data,2), 'ax', ax,...
    'mask', params.wm_mask,'clip_vals', false, 'limits',diff_limits);
set(ax,'CameraViewAngle',6.028);
set(ax,'CameraUpVector',subj_info.camera_up_vector(view));
set(ax,'CameraPosition',subj_info.camera_position(view));
freezeColors;

%ax=subplot(2,4,7);
ax=subplot(1,4,3);
%ax=0;
[ax,~]=plot_surface_metric(pial_inflated_mesh, mean(pial_wm_diff_data,2), 'ax', ax,...
    'mask', params.mask,'clip_vals', false, 'limits',pial_wm_diff_limits);
set(ax,'CameraViewAngle',6.028);
set(ax,'CameraUpVector',subj_info.camera_up_vector(view));
set(ax,'CameraPosition',subj_info.camera_position(view));
freezeColors;

% ax=subplot(2,4,8);
% %figure();
% nbins=100.0;
% masked_pial_wm_diff=pial_wm_diff_data;
% if length(params.mask)
%     masked_pial_wm_diff=pial_wm_diff_data(params.mask);
% end
% bin_width=(max(masked_pial_wm_diff)-min(masked_pial_wm_diff))/nbins;
% [f,x]=hist(masked_pial_wm_diff,[min(masked_pial_wm_diff):bin_width:max(masked_pial_wm_diff)]);      
% bar(x,f/sum(f).*100.0,'b');
% xlim(pial_wm_diff_limits);
% %if params.hist_ylim>0
% %    ylim([0 params.hist_ylim]);
% %end
% xlabel('Pial-White');
% ylabel('Percent vertices');

if length(params.output_file)>0
    saveas(fig, params.output_file, params.output_format);
end
