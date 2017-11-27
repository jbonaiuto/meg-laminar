function compare_session_layers(subj_info, session_num, contrast, foi_dir, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'save_results', true,...
    'inv_type','EBB','patch_size', 0.4, 'surf_dir','');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end
if length(params.surf_dir)==0
    params.surf_dir=fullfile(params.data_dir,'surf');
end

% Get map from white matter to pial surface
orig_white_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),...
    'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),...
    'surf','ds_white.hires.deformed.surf.gii');
orig_pial_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),...
    'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),...
    'surf','ds_pial.hires.deformed.surf.gii');

pial_white_map=map_pial_to_white(white_mesh, pial_mesh, 'mapType', 'link',...
    'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh);

nvertices=length(pial_white_map);

baseline_dir=fullfile(foi_dir, ...
    sprintf('t%d_%d', contrast.baseline_woi(1), contrast.baseline_woi(2)));
pial_baseline=load_woi_trials(baseline_dir,...
    sprintf('pial_br%s_%d_1_', subj_info.subj_id, session_num),...
    contrast.baseline_woi, contrast.foi, nvertices);
white_baseline=load_woi_trials(baseline_dir,...
    sprintf('white_br%s_%d_1_', subj_info.subj_id, session_num),...
    contrast.baseline_woi, contrast.foi, nvertices);
    
woi_dir=fullfile(foi_dir, sprintf('t%d_%d', contrast.comparison_woi(1),...
    contrast.comparison_woi(2)));
pial_vals=load_woi_trials(woi_dir, sprintf('pial_br%s_%d_1_',...
    subj_info.subj_id, session_num), contrast.comparison_woi,...
    contrast.foi, nvertices);
white_vals=load_woi_trials(woi_dir, sprintf('white_br%s_%d_1_',...
    subj_info.subj_id, session_num), contrast.comparison_woi,...
    contrast.foi, nvertices);

pial_diff=pial_vals-pial_baseline;
white_diff=white_vals-white_baseline;

% Save pial/wm diff
write_metric_gifti(fullfile(foi_dir, sprintf('pial.%s.diff.dat',...
    contrast.comparison_name)), pial_diff);
write_metric_gifti(fullfile(foi_dir, sprintf('white.%s.diff.dat',...
    contrast.comparison_name)), white_diff);

% Compare pial values at two wois
[H,pvals,ci,STATS]=ttest(pial_diff');
pial_tvals=STATS.tstat';

% Save pial comparison
write_metric_gifti(fullfile(foi_dir, sprintf('pial.%s.t.dat',...
    contrast.comparison_name)), pial_tvals);

pial_d=(mean(pial_vals,2)-mean(pial_baseline,2))./sqrt((var(pial_vals,[],2)+var(pial_baseline,[],2))/2);
write_metric_gifti(fullfile(foi_dir, sprintf('pial.%s.d.dat',...
    contrast.comparison_name)), pial_d);

% Compare white matter values at two wois
[H,pvals,ci,STATS]=ttest(white_diff');
white_tvals=STATS.tstat';

% Save white matter comparison
write_metric_gifti(fullfile(foi_dir, sprintf('white.%s.t.dat',...
    contrast.comparison_name)), white_tvals);        

white_d=(mean(white_vals,2)-mean(white_baseline,2))./sqrt((var(white_vals,[],2)+var(white_baseline,[],2))/2);
write_metric_gifti(fullfile(foi_dir, sprintf('white.%s.d.dat',...
    contrast.comparison_name)), white_d);

% Compare pial and white matter differences
pial_white_diff=abs(pial_diff)-abs(white_diff(pial_white_map,:));
[H,pvals,ci,STATS]=ttest(pial_white_diff');
pial_white_tvals=STATS.tstat';

% Save pial - white matter diff
write_metric_gifti(fullfile(foi_dir, sprintf('pial-white.%s.diff.dat',...
    contrast.comparison_name)), pial_white_diff);        
write_metric_gifti(fullfile(foi_dir, sprintf('pial-white.%s.t.dat',...
    contrast.comparison_name)), pial_white_tvals);        

