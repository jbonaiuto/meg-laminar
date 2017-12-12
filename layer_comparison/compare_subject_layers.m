function compare_subject_layers(subj_info, contrast, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'save_results', true,...
    'inv_type','EBB','patch_size', 5.0, 'surf_dir','', 'filter_sessions', true);  %define default values
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
orig_white_mesh=fullfile(params.surf_dir,subj_info.subj_id,'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,subj_info.subj_id,'surf','ds_white.hires.deformed.surf.gii');
orig_pial_mesh=fullfile(params.surf_dir,subj_info.subj_id,'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,subj_info.subj_id,'surf','ds_pial.hires.deformed.surf.gii');

pial_white_map=map_pial_to_white(white_mesh, pial_mesh, 'mapType', 'link',...
    'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh);

nvertices=length(pial_white_map);

pial_baseline=[];
white_baseline=[];

for session_num=1:length(subj_info.sessions)
    foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id, num2str(session_num), 'grey_coreg', params.inv_type, ['p' num2str(params.patch_size)], contrast.zero_event, ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);

    baseline_dir=fullfile(foi_dir,['t' num2str(contrast.baseline_woi(1)) '_' num2str(contrast.baseline_woi(2))]);
    session_pial_baseline=load_woi_trials(baseline_dir, ['pial_br' subj_info.subj_id '_' num2str(session_num) '_1_'], contrast.baseline_woi, contrast.foi, nvertices);
    session_white_baseline=load_woi_trials(baseline_dir, ['white_br' subj_info.subj_id '_' num2str(session_num) '_1_'], contrast.baseline_woi, contrast.foi, nvertices);
    ntrials=size(session_pial_baseline,2);
    pial_baseline(:,end+1:end+ntrials)=session_pial_baseline;
    white_baseline(:,end+1:end+ntrials)=session_white_baseline;
end

pial_vals=[];
white_vals=[];
for session_num=1:length(subj_info.sessions)
    foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id, num2str(session_num), 'grey_coreg', params.inv_type, ['p' num2str(params.patch_size)], contrast.zero_event, ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
    woi_dir=fullfile(foi_dir,['t' num2str(contrast.comparison_woi(1)) '_' num2str(contrast.comparison_woi(2))]);
    session_pial_vals=load_woi_trials(woi_dir, ['pial_br' subj_info.subj_id '_' num2str(session_num) '_1_'], contrast.comparison_woi, contrast.foi, nvertices);
    session_white_vals=load_woi_trials(woi_dir, ['white_br' subj_info.subj_id '_' num2str(session_num) '_1_'], contrast.comparison_woi, contrast.foi, nvertices);
    ntrials=size(session_pial_vals,2);
    pial_vals(:,end+1:end+ntrials)=session_pial_vals;
    white_vals(:,end+1:end+ntrials)=session_white_vals;
end

foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id, 'grey_coreg', params.inv_type, ['p' num2str(params.patch_size)], contrast.zero_event, ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
mkdir(foi_dir);

pial_diff=pial_vals-pial_baseline;
white_diff=white_vals-white_baseline;

% Save pial/wm diff
write_metric_gifti(fullfile(foi_dir, ['pial.' contrast.comparison_name '.diff.dat']), pial_diff);
write_metric_gifti(fullfile(foi_dir, ['white.' contrast.comparison_name '.diff.dat']), white_diff);

% Compare pial values at two wois
[H,pvals,ci,STATS]=ttest(pial_diff');
pial_tvals=STATS.tstat';

% Save pial comparison
write_metric_gifti(fullfile(foi_dir, ['pial.' contrast.comparison_name '.t.dat']), pial_tvals);

pial_d=(mean(pial_vals,2)-mean(pial_baseline,2))./sqrt((var(pial_vals,[],2)+var(pial_baseline,[],2))/2);
write_metric_gifti(fullfile(foi_dir, ['pial.' contrast.comparison_name '.d.dat']), pial_d);

% Compare white matter values at two wois
[H,pvals,ci,STATS]=ttest(white_diff');
white_tvals=STATS.tstat';

% Save white matter comparison
write_metric_gifti(fullfile(foi_dir, ['white.' contrast.comparison_name '.t.dat']), white_tvals);        

white_d=(mean(white_vals,2)-mean(white_baseline,2))./sqrt((var(white_vals,[],2)+var(white_baseline,[],2))/2);
write_metric_gifti(fullfile(foi_dir, ['white.' contrast.comparison_name '.d.dat']), white_d);

% Compare pial and white matter differences
pial_white_diff=abs(pial_diff)-abs(white_diff(pial_white_map,:));
[H,pvals,ci,STATS]=ttest(pial_white_diff');
pial_white_tvals=STATS.tstat';

% Save pial - white matter diff
write_metric_gifti(fullfile(foi_dir, ['pial-white.' contrast.comparison_name '.diff.dat']), pial_white_diff);        
write_metric_gifti(fullfile(foi_dir, ['pial-white.' contrast.comparison_name '.t.dat']), pial_white_tvals);        


