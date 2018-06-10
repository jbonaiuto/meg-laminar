function [pial_mask,wm_mask,mask]=compute_roi(subj_info, foi_dir, comparison_name, ...
    thresh_type, pial_name, wm_name, pial_inflated_name, wm_inflated_name,...
    pial_white_map, white_pial_map, lfn_filename, varargin)

% Parse inputs
defaults = struct('hemisphere','','pial_hemisphere_map',[],...
    'white_hemisphere_map',[],'region','','thresh_percentile',95,...
    'type','mean','origPial','','origWhite','','recompute',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

pial=gifti(pial_name);
wm=gifti(wm_name);

nvertices=size(pial.vertices,1);
    
D=spm_eeg_load(lfn_filename);
pial_lfn=sqrt(sum(D.inv{1}.inverse.L(:,nvertices+1:end).^2,1))';
pial_lfn=spm_mesh_smooth(pial, pial_lfn, 20);
wm_lfn=sqrt(sum(D.inv{1}.inverse.L(:,1:nvertices).^2,1))';
wm_lfn=spm_mesh_smooth(wm, wm_lfn, 20);        
mapped_wm_lfn=wm_lfn(pial_white_map);

switch params.type
    case 't'
        pial_t=gifti(fullfile(foi_dir,['pial.' comparison_name '.t.gii']));
        pial_metric=pial_t.cdata(:);
        wm_t=gifti(fullfile(foi_dir,['white.' comparison_name '.t.gii']));
        wm_metric=wm_t.cdata(:);        
    case 'd'
        pial_d=gifti(fullfile(foi_dir,['pial.' comparison_name '.d.gii']));
        pial_metric=pial_d.cdata(:);
        wm_d=gifti(fullfile(foi_dir,['white.' comparison_name '.d.gii']));
        wm_metric=wm_d.cdata(:);  
    case 'mean'
        pial_mean=gifti(fullfile(foi_dir,['pial.' comparison_name '.diff.gii']));
        pial_metric=mean(pial_mean.cdata(:,:),2);
        wm_mean=gifti(fullfile(foi_dir,['white.' comparison_name '.diff.gii']));
        wm_metric=mean(wm_mean.cdata(:,:),2);  
end
mapped_wm_metric=wm_metric(pial_white_map);

prethresh_pial_mask=intersect(find(~isinf(pial_metric)), find(pial_lfn>22));
prethresh_wm_mask=intersect(find(~isinf(wm_metric)), find(wm_lfn>19));
prethresh_mapped_wm_mask=intersect(find(~isinf(mapped_wm_metric)), find(mapped_wm_lfn>19));

switch params.hemisphere
    case 'left'
        prethresh_pial_mask=intersect(prethresh_pial_mask, find(params.pial_hemisphere_map==1));
        prethresh_wm_mask=intersect(prethresh_wm_mask, find(params.white_hemisphere_map==1));
        prethresh_mapped_wm_mask=intersect(prethresh_mapped_wm_mask, find(params.white_hemisphere_map(pial_white_map)==1));
    case 'right'
        prethresh_pial_mask=intersect(prethresh_pial_mask, find(params.pial_hemisphere_map==2));
        prethresh_wm_mask=intersect(prethresh_wm_mask, find(params.white_hemisphere_map==2));
        prethresh_mapped_wm_mask=intersect(prethresh_mapped_wm_mask, find(params.white_hemisphere_map(pial_white_map)==2));
end
    
if length(params.region)>0
    switch params.region
        case 'visual'
            [lh_pial_mask,lh_wm_mask]=get_anatomical_mask('lh_visual', pial_name, pial_inflated_name,...
                wm_name, wm_inflated_name, subj_info.coords('lh_visual'), 50,...
                'origPial', params.origPial, 'origWhite', params.origWhite, 'recompute', params.recompute);
            [rh_pial_mask,rh_wm_mask]=get_anatomical_mask('rh_visual', pial_name, pial_inflated_name,...
                wm_name, wm_inflated_name, subj_info.coords('rh_visual'), 50,...
                'origPial', params.origPial, 'origWhite', params.origWhite, 'recompute', params.recompute);
            prethresh_pial_mask=intersect(prethresh_pial_mask, union(lh_pial_mask,rh_pial_mask));
            prethresh_wm_mask=intersect(prethresh_wm_mask, union(lh_wm_mask,rh_wm_mask));
            prethresh_mapped_wm_mask=intersect(prethresh_mapped_wm_mask, white_pial_map(prethresh_wm_mask));
        otherwise
            [anat_pial_mask,anat_wm_mask]=get_anatomical_mask(params.region, pial_name, pial_inflated_name,...
                wm_name, wm_inflated_name, subj_info.coords(params.region), 75,...
                'origPial', params.origPial, 'origWhite', params.origWhite, 'recompute', params.recompute);
            prethresh_pial_mask=intersect(prethresh_pial_mask, anat_pial_mask);
            prethresh_wm_mask=intersect(prethresh_wm_mask, anat_wm_mask);
            prethresh_mapped_wm_mask=intersect(prethresh_mapped_wm_mask, white_pial_map(anat_wm_mask));
    end
end

% Find vertices greater than threshold
if strcmp(thresh_type,'lower')
    pial_threshold=prctile(pial_metric(prethresh_pial_mask),params.thresh_percentile);
    wm_threshold=prctile(wm_metric(prethresh_wm_mask),params.thresh_percentile);
    % Create pial and white masks and mapped white mask
    pial_mask=intersect(prethresh_pial_mask, find(pial_metric>pial_threshold));
    wm_mask=intersect(prethresh_wm_mask, find(wm_metric>wm_threshold));
    mapped_wm_mask=intersect(prethresh_mapped_wm_mask, find(mapped_wm_metric>wm_threshold));
% Find vertices less than threshold
else
    pial_threshold=prctile(pial_metric(prethresh_pial_mask),100-params.thresh_percentile);
    wm_threshold=prctile(wm_metric(prethresh_wm_mask),100-params.thresh_percentile);
    % Create pial and white maskss and mapped white mask
    pial_mask=intersect(prethresh_pial_mask, find(pial_metric<pial_threshold));
    wm_mask=intersect(prethresh_wm_mask, find(wm_metric<wm_threshold));        
    mapped_wm_mask=intersect(prethresh_mapped_wm_mask, find(mapped_wm_metric<wm_threshold));        
end
mask=union(pial_mask, mapped_wm_mask);