function split_inversion_results(subj_info, grey_coreg_dir, foi, woi, varargin)

% Parse inputs
defaults = struct('patch_size',0.4);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

woi_dir=fullfile(grey_coreg_dir, ['p' num2str(params.patch_size)], ['f' num2str(foi(1)) '_' num2str(foi(2))], ['t' num2str(woi(1)) '_' num2str(woi(2))]);

[files,dirs] = spm_select('List', woi_dir, ['^r' subj_info.subj_id '.*\.gii']);
for f=1:size(files,1)
    filename=deblank(files(f,:));
    parts=splitstring(filename,'.');
    prefix=parts{1};
    trial_mesh=gifti(fullfile(woi_dir,filename));
    n_combined_vertices=length(trial_mesh.cdata);
    n_vertices=round(n_combined_vertices/2);

    c=file_array(fullfile(woi_dir,['pial_' prefix '.dat']),[n_vertices 1],'FLOAT32-LE',0,1,0);
    c(:)=trial_mesh.cdata(n_vertices+1:end);
    pial_surf = gifti;
    pial_surf.cdata = c;
    save(pial_surf, fullfile(woi_dir,['pial_' prefix '.gii']), 'ExternalFileBinary');

    c=file_array(fullfile(woi_dir,['white_' prefix '.dat']),[n_vertices 1],'FLOAT32-LE',0,1,0);
    c(:)=trial_mesh.cdata(1:n_vertices);
    white_surf = gifti;
    white_surf.cdata = c;
    save(white_surf, fullfile(woi_dir,['white_' prefix '.gii']), 'ExternalFileBinary');

    delete(fullfile(woi_dir,[prefix '.dat']),fullfile(woi_dir,[prefix '.gii']));
end
