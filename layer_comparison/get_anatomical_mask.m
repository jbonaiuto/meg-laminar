function [pial_mask,wm_mask]=get_anatomical_mask(roi_name, pial_name, ...
    inflated_pial_name, white_name, inflated_white_name, inflated_coord,...
    inflated_radius, varargin)
% 1=left, 2=right

% Parse inputs
defaults = struct('recompute',false,'origPial','','origWhite','');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

[meshpath pial_mesh_file ext]=fileparts(pial_name);
pial_mask_file=fullfile(meshpath, sprintf('%s_%s.mat', roi_name, pial_mesh_file));
[meshpath wm_mesh_file ext]=fileparts(white_name);
wm_mask_file=fullfile(meshpath, sprintf('%s_%s.mat', roi_name, wm_mesh_file));
if exist(pial_mask_file,'file')~=2 || exist(wm_mask_file,'file')~=2 || params.recompute
    disp('Recomputing mask');
    
    pial_white_map=map_pial_to_white(white_name, pial_name, 'mapType', 'link',...
        'origPial', params.origPial, 'origWhite', params.origWhite);

    pial=gifti(pial_name);
    n_pial_vertices=size(pial.vertices,1);
    inflated_pial_surface=gifti(inflated_pial_name);
    pial_dist=sqrt(sum((inflated_pial_surface.vertices-repmat(inflated_coord,n_pial_vertices,1)).^2,2));
    pial_mask=find(pial_dist<=inflated_radius);
    save(pial_mask_file,'pial_mask');
    
    white=gifti(white_name);
    n_wm_vertices=size(white.vertices,1);
    inflated_white_surface=gifti(inflated_white_name);
    wm_coord=inflated_white_surface.vertices(pial_white_map(find(pial_dist==min(pial_dist))),:);
    wm_dist=sqrt(sum((inflated_white_surface.vertices-repmat(wm_coord,n_wm_vertices,1)).^2,2));
    wm_mask=find(wm_dist<=inflated_radius);
    save(wm_mask_file,'wm_mask');
else
    a=load(pial_mask_file);
    pial_mask=a.pial_mask;
    a=load(wm_mask_file);
    wm_mask=a.wm_mask;
end
