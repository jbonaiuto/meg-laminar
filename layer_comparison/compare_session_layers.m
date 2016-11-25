function compare_session_layers(subj_info, session_num, foi, woi, baseline, comparison_name, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'save_results', true, 'patch_size', 0.4, 'white_pial_map', [],'surf_dir','');  %define default values
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
if length(params.white_pial_map)==0
    white=gifti(fullfile(params.surf_dir, [subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii'));
    pial=gifti(fullfile(params.surf_dir, [subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii'));
    params.white_pial_map=map_white_to_pial(white, pial);
end

foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id, num2str(session_num), 'grey_coreg', ['p' num2str(params.patch_size)], ['f' num2str(foi(1)) '_' num2str(foi(2))]);

load(fullfile(foi_dir,sprintf('r%s_%d.mat', subj_info.subj_id, session_num)));
ntrials=length(D.trials);
    
spm('defaults', 'EEG');

woi_dir=fullfile(foi_dir, ['t' num2str(woi(1)) '_' num2str(woi(2))]);
baseline_dir=fullfile(foi_dir, ['t' num2str(baseline(1)) '_' num2str(baseline(2))]);

% Load all pial data from woi
pial_woi_trials=[];
%[files,dirs] = spm_select('List', woi_dir, ['^pial_r' subj_info.subj_id '_' num2str(session_num) '_1_t' num2str(woi(1)) '_' num2str(woi(2)) '_f' num2str(foi(1)) '_' num2str(foi(2)) '.*\.gii']);
for t=1:ntrials
    filename=fullfile(woi_dir, sprintf('pial_br%s_%d_1_t%d_%d_f%d_%d_1_%d.gii', subj_info.subj_id, session_num, woi(1), woi(2), foi(1), foi(2), t));
    trial_mesh=gifti(filename);
    pial_woi_trials(:,t)=trial_mesh.cdata(:);
end

% Load all pial data from baseline
pial_baseline_trials=[];
for t=1:ntrials
    filename=fullfile(baseline_dir, sprintf('pial_br%s_%d_1_t%d_%d_f%d_%d_1_%d.gii', subj_info.subj_id, session_num, baseline(1), baseline(2), foi(1), foi(2), t));
    trial_mesh=gifti(filename);
    pial_baseline_trials(:,t)=trial_mesh.cdata(:);
end
% pial_baseline=repmat(mean(pial_baseline_trials,2),1,ntrials);
% pial_diff=(pial_woi_trials-pial_baseline)./pial_baseline;
%pial_diff=pial_woi_trials;
pial_diff=pial_woi_trials-pial_baseline_trials;

% Load all white matter data from woi
white_woi_trials=[];
%[files,dirs] = spm_select('List', woi_dir, ['^white_r' subj_info.subj_id '_' num2str(session_num) '_1_t' num2str(woi(1)) '_' num2str(woi(2)) '_f' num2str(foi(1)) '_' num2str(foi(2)) '.*\.gii']);
for t=1:ntrials
    filename=fullfile(woi_dir, sprintf('white_br%s_%d_1_t%d_%d_f%d_%d_1_%d.gii', subj_info.subj_id, session_num, woi(1), woi(2), foi(1), foi(2), t));
    trial_mesh=gifti(filename);
    white_woi_trials(:,t)=trial_mesh.cdata(:);
end

% Load all white matter data from baseline
white_baseline_trials=[];
%[files,dirs] = spm_select('List', baseline_dir, ['^white_r' subj_info.subj_id '_' num2str(session_num) '_1_t' num2str(baseline(1)) '_' num2str(baseline(2)) '_f' num2str(foi(1)) '_' num2str(foi(2)) '.*\.gii']);
for t=1:ntrials
    filename=fullfile(baseline_dir, sprintf('white_br%s_%d_1_t%d_%d_f%d_%d_1_%d.gii', subj_info.subj_id, session_num, baseline(1), baseline(2), foi(1), foi(2), t));
    trial_mesh=gifti(filename);
    white_baseline_trials(:,t)=trial_mesh.cdata(:);
end
% white_baseline=repmat(mean(white_baseline_trials,2),1,ntrials);
% white_diff=(white_woi_trials-white_baseline)./white_baseline;
%white_diff=white_woi_trials;
white_diff=white_woi_trials-white_baseline_trials;

if params.save_results
    % Save pial diff
    c=file_array(fullfile(foi_dir, ['pial.' comparison_name '.diff.shape.dat']),size(pial_diff,1),'FLOAT32-LE',0,1,0);
    c(:)=mean(pial_diff,2);
    pial_diff_surf = gifti;
    pial_diff_surf.cdata = c;
    save(pial_diff_surf, fullfile(foi_dir, ['pial.' comparison_name '.diff.shape.gii']), 'ExternalFileBinary');

    % Compare pial values at two wois
    [H,pvals,ci,STATS]=ttest(pial_diff');
    pial_tvals=STATS.tstat';

    % Save pial comparison
    c=file_array(fullfile(foi_dir, ['pial.' comparison_name '.t.shape.dat']),size(pial_tvals),'FLOAT32-LE',0,1,0);
    c(:)=pial_tvals;
    pial_t_surf = gifti;
    pial_t_surf.cdata = c;
    save(pial_t_surf, fullfile(foi_dir, ['pial.' comparison_name '.t.shape.gii']), 'ExternalFileBinary');

    % Save white matter diff
    c=file_array(fullfile(foi_dir, ['white.' comparison_name '.diff.shape.dat']),size(white_diff,1),'FLOAT32-LE',0,1,0);
    c(:)=mean(white_diff,2);
    white_diff_surf = gifti;
    white_diff_surf.cdata = c;
    save(white_diff_surf, fullfile(foi_dir, ['white.' comparison_name '.diff.shape.gii']), 'ExternalFileBinary');

    % Compare white matter values at two wois
    [H,pvals,ci,STATS]=ttest(white_diff');
    white_tvals=STATS.tstat';

    % Save white matter comparison
    c=file_array(fullfile(foi_dir, ['white.' comparison_name '.t.shape.dat']),size(white_tvals),'FLOAT32-LE',0,1,0);
    c(:)=white_tvals;
    white_t_surf = gifti;
    white_t_surf.cdata = c;
    save(white_t_surf, fullfile(foi_dir, ['white.' comparison_name '.t.shape.gii']), 'ExternalFileBinary');

    % Compare pial and white matter differences
    pial_white_diff=abs(pial_diff)-abs(white_diff(params.white_pial_map,:));
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
