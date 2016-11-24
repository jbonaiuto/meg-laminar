function [pial_perc_change, white_perc_change]=compare_session_layers(subj_info, session_num, foi, woi, baseline, comparison_name, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'save_results', true, 'patch_size', 0.4, 'white_pial_map', []);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Get map from white matter to pial surface
if length(params.white_pial_map)==0
    white=gifti(fullfile(params.data_dir, 'surf',[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii'));
    pial=gifti(fullfile(params.data_dir, 'surf',[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii'));
    params.white_pial_map=map_white_to_pial(white, pial);
end

foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id, num2str(session_num), 'grey_coreg', ['p' num2str(params.patch_size)], ['f' num2str(foi(1)) '_' num2str(foi(2))]);

spm('defaults', 'EEG');

woi_dir=fullfile(foi_dir, ['t' num2str(woi(1)) '_' num2str(woi(2))]);
baseline_dir=fullfile(foi_dir, ['t' num2str(baseline(1)) '_' num2str(baseline(2))]);

% Load all pial data from woi
pial_woi_trials=[];
[files,dirs] = spm_select('List', woi_dir, ['^pial_r' subj_info.subj_id '_' num2str(session_num) '_1_t' num2str(woi(1)) '_' num2str(woi(2)) '_f' num2str(foi(1)) '_' num2str(foi(2)) '.*\.gii']);
for f=1:size(files,1)
    filename=deblank(files(f,:));
    trial_mesh=gifti(fullfile(woi_dir,filename));
    pial_woi_trials(:,end+1)=trial_mesh.cdata(:);
end

% Load all pial data from baseline
pial_baseline_trials=[];
[files,dirs] = spm_select('List', baseline_dir, ['^pial_r' subj_info.subj_id '_' num2str(session_num) '_1_t' num2str(baseline(1)) '_' num2str(baseline(2)) '_f' num2str(foi(1)) '_' num2str(foi(2)) '.*\.gii']);
for f=1:size(files,1)
    filename=deblank(files(f,:));
    trial_mesh=gifti(fullfile(baseline_dir,filename));
    pial_baseline_trials(:,end+1)=trial_mesh.cdata(:);
end
% Compute percentage change from baseline in pial surface
% TODO - this is per trial baseline, what about mean baseline?
pial_perc_change=(pial_woi_trials-pial_baseline_trials)./pial_baseline_trials;

% Load all white matter data from woi
white_woi_trials=[];
[files,dirs] = spm_select('List', woi_dir, ['^white_r' subj_info.subj_id '_' num2str(session_num) '_1_t' num2str(woi(1)) '_' num2str(woi(2)) '_f' num2str(foi(1)) '_' num2str(foi(2)) '.*\.gii']);
for f=1:size(files,1)
    filename=deblank(files(f,:));
    trial_mesh=gifti(fullfile(woi_dir,filename));
    white_woi_trials(:,end+1)=trial_mesh.cdata(:);
end

% Load all white matter data from baseline
white_baseline_trials=[];
[files,dirs] = spm_select('List', baseline_dir, ['^white_r' subj_info.subj_id '_' num2str(session_num) '_1_t' num2str(baseline(1)) '_' num2str(baseline(2)) '_f' num2str(foi(1)) '_' num2str(foi(2)) '.*\.gii']);
for f=1:size(files,1)
    filename=deblank(files(f,:));
    trial_mesh=gifti(fullfile(baseline_dir, filename));
    white_baseline_trials(:,end+1)=trial_mesh.cdata(:);
end
% Compute percent change from baseline in white matter surface
% TODO - this is per trial baseline, what about mean baseline?
white_perc_change=(white_woi_trials-white_baseline_trials)./white_baseline_trials;

if params.save_results
    % Save pial percent change
    c=file_array(fullfile(foi_dir, ['pial.' comparison_name '.perc_change.shape.dat']),size(pial_perc_change,1),'FLOAT32-LE',0,1,0);
    c(:)=mean(pial_perc_change,2);
    pial_perc_change_surf = gifti;
    pial_perc_change_surf.cdata = c;
    save(pial_perc_change_surf, fullfile(foi_dir, ['pial.' comparison_name '.perc_change.shape.gii']), 'ExternalFileBinary');

    % Compute significance of change from baseline in pial
    [H,pvals,ci,STATS]=ttest(pial_perc_change');
    pial_tvals=STATS.tstat';

    % Save pial comparison
    c=file_array(fullfile(foi_dir, ['pial.' comparison_name '.t.shape.dat']),size(pial_tvals),'FLOAT32-LE',0,1,0);
    c(:)=pial_tvals;
    pial_t_surf = gifti;
    pial_t_surf.cdata = c;
    save(pial_t_surf, fullfile(foi_dir, ['pial.' comparison_name '.t.shape.gii']), 'ExternalFileBinary');

    % Save white matter percent change
    c=file_array(fullfile(foi_dir, ['white.' comparison_name '.perc_change.shape.dat']),size(white_perc_change,1),'FLOAT32-LE',0,1,0);
    c(:)=mean(white_perc_change,2);
    white_perc_change_surf = gifti;
    white_perc_change_surf.cdata = c;
    save(white_perc_change_surf, fullfile(foi_dir, ['white.' comparison_name '.perc_change.shape.gii']), 'ExternalFileBinary');

    % Compute significance of change from baseline in white matter
    [H,pvals,ci,STATS]=ttest(white_perc_change');
    white_tvals=STATS.tstat';

    % Save white matter comparison
    c=file_array(fullfile(foi_dir, ['white.' comparison_name '.t.shape.dat']),size(white_tvals),'FLOAT32-LE',0,1,0);
    c(:)=white_tvals;
    white_t_surf = gifti;
    white_t_surf.cdata = c;
    save(white_t_surf, fullfile(foi_dir, ['white.' comparison_name '.t.shape.gii']), 'ExternalFileBinary');

    % Compare pial and white matter differences
    pial_white_diff=abs(pial_perc_change)-abs(white_perc_change(params.white_pial_map,:));
    [H,pvals,ci,STATS]=ttest(pial_white_diff');
    pial_white_tvals=STATS.tstat';

    % Save pial - white matter diff
    c=file_array(fullfile(foi_dir, ['pial-white.' comparison_name '.diff.shape.dat']),size(pial_white_diff,1),'FLOAT32-LE',0,1,0);
    c(:)=mean(pial_white_diff,2);
    pial_white_diff_surf = gifti;
    pial_white_diff_surf.cdata = c;
    save(pial_white_diff_surf, fullfile(foi_dir,['pial-white.' comparison_name '.diff.shape.gii']), 'ExternalFileBinary');

    % Save pial - white matter comparison
    c=file_array(fullfile(foi_dir, ['pial-white.' comparison_name '.t.shape.dat']),size(pial_white_tvals),'FLOAT32-LE',0,1,0);
    c(:)=pial_white_tvals;
    pial_white_t_surf = gifti;
    pial_white_t_surf.cdata = c;
    save(pial_white_t_surf, fullfile(foi_dir,['pial-white.' comparison_name '.t.shape.gii']), 'ExternalFileBinary');

end
