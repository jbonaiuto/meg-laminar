function hemisphere_map=get_hemisphere_map(surface_name, orig_surface_name, varargin)
% 1=left, 2=right

% Parse inputs
defaults = struct('recompute',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

[meshpath mesh_file ext]=fileparts(surface_name);
hemi_file=fullfile(meshpath, sprintf('hemisphere_map_%s.mat', mesh_file));
if exist(hemi_file,'file')~=2 || params.recompute
    disp('Recomputing hemisphere map');
    surface=gifti(surface_name);
    n_vertices=size(surface.vertices,1);
    hemisphere_map=ones(1,n_vertices);
    orig_surface=gifti(orig_surface_name);
    %n_orig_vertices=size(orig_surface.vertices,1);
    [meshpath orig_mesh_file ext]=fileparts(orig_surface_name);
    lh_orig_surface=gifti(fullfile(meshpath,sprintf('lh.%s.gii',orig_mesh_file)));
    n_lh_vertices=size(lh_orig_surface.vertices,1);
    % Map from downsampled pial to original pial surface
    surf_ds_surf_map=dsearchn(orig_surface.vertices,surface.vertices);
    hemisphere_map(find(surf_ds_surf_map>n_lh_vertices))=2;
    save(hemi_file,'hemisphere_map');
else
    a=load(hemi_file);
    hemisphere_map=a.hemisphere_map;
end
