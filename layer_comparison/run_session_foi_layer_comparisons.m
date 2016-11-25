% Comparison wois - matrix, nx2, n=number of comparisons
% Baseline wois - matrix, nx2, n=number of comparisons
% Resp-locked Beta ERD, rebound:
%     zero_event='resp'
%     foi=[15 30]
%     inversion_woi=[-2000 2000]
%     comparison_names={'resp_beta_erd','resp_beta_rebound'}
%     comparison_wois=[-250 250;500 1000];
%     baseline_wois=[-1500 -1000;-1500 -1000]
% Resp-locked gamma MRGS
%     zero_event='resp'
%     foi=[60 90]
%     inversion_woi=[-1500 1500]
%     comparison_names={'resp_mrgs'}
%     comparison_wois=[-100 200];
%     baseline_wois=[-1000 -700]
function run_session_foi_layer_comparisons(subj_info, session_num, zero_event, foi, inversion_woi, comparison_names, comparison_wois, baseline_woi, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'patch_size',0.4, 'surf_dir', '', 'mri_dir', '');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end
if length(params.surf_dir)==0
    params.surf_dir=fullfile(params.data_dir,'surf');
end

invert_grey(subj_info, session_num, zero_event, foi, inversion_woi, baseline_woi, 'data_dir', params.data_dir, 'patch_size', params.patch_size, 'surf_dir', params.surf_dir, 'mri_dir', params.mri_dir);

all_wois=[comparison_wois; baseline_woi];
unique_wois=unique(all_wois,'rows');

for i=1:size(unique_wois,1)
    extract_inversion_source(subj_info, session_num, foi, unique_wois(i,:), 'data_dir', params.data_dir, 'patch_size', params.patch_size);
end

white=gifti(fullfile(params.surf_dir, [subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii'));
pial=gifti(fullfile(params.surf_dir, [subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii'));
white_pial_map=map_white_to_pial(white, pial);

for i=1:length(comparison_names)
    comparison_name=comparison_names{i};
    compare_session_layers(subj_info, session_num, foi, comparison_wois(i,:), baseline_woi, comparison_name, 'data_dir', params.data_dir, 'patch_size', params.patch_size, 'white_pial_map', white_pial_map, 'surf_dir', params.surf_dir);
end
