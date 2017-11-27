function run_all_subject_coregerr(subj_info, contrasts, varargin)

defaults = struct('data_dir', 'd:/pred_coding', 'inv_type', 'EBB',...
    'patch_size',0.4, 'surf_dir', '', 'filter_sessions',true, 'iterations',10,...
    'shift_magnitude', 10);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end
if length(params.surf_dir)==0
    params.surf_dir=fullfile(params.data_dir,'surf');
end

for idx=1:params.iterations
    for i=1:length(contrasts)
        contrast=contrasts(i);
        compare_subject_layers_coregerr(subj_info, contrast, idx,...
            'data_dir', params.data_dir, 'save_results', true, 'inv_type', params.inv_type, ...
            'patch_size', params.patch_size, 'surf_dir', params.surf_dir,...
            'filter_sessions', params.filter_sessions, 'shift_magnitude', params.shift_magnitude);
    end        
end
