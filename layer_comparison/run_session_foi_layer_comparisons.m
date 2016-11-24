% Comparison wois - matrix, nx2, n=number of comparisons
% Baseline wois - matrix, nx2, n=number of comparisons
function run_session_foi_layer_comparisons(subj_info, session_num, zero_event, foi, inversion_woi, comparison_names, comparison_wois, baseline_wois, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'patch_size',0.4);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

invert_grey(subj_info, session_num, zero_event, foi, inversion_woi, 'data_dir', params.data_dir, 'patch_size', params.patch_size);

all_wois=[comparison_wois; baseline_wois];
unique_wois=unique(all_wois,'rows');

extract_inversion_source(subj_info, session_num, foi, unique_wois, 'data_dir', params.data_dir, 'patch_size', params.patch_size);

white=gifti(fullfile(params.data_dir, 'surf',[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii'));
pial=gifti(fullfile(params.data_dir, 'surf',[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii'));
white_pial_map=map_white_to_pial(white, pial);

for i=1:length(comparison_names)
    comparison_name=comparison_name{i};
    compare_session_layers(subj_info, session_num, foi, comparison_wois(i,:), baseline_wois(i,:), comparison_name, 'data_dir', params.data_dir, 'patch_size', params.patch_size, 'white_pial_map', white_pial_map);
end
