function [pial_diff, white_diff]=compare_session_layers(subj_info, session_num, foi, woi1, woi2, comparison_name, varargin)

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

woi1_dir=fullfile(foi_dir, ['t' num2str(woi1(1)) '_' num2str(woi1(2))]);
woi2_dir=fullfile(foi_dir, ['t' num2str(woi2(1)) '_' num2str(woi2(2))]);

% Load all pial data from woi1
pial_woi1_trials=[];
[files,dirs] = spm_select('List', woi1_dir, ['^pial_r' subj_info.subj_id '_' num2str(session_num) '_1_t' num2str(woi1(1)) '_' num2str(woi1(2)) '_f' num2str(foi(1)) '_' num2str(foi(2)) '.*\.gii']);
for f=1:size(files,1)
    filename=deblank(files(f,:));
    trial_mesh=gifti(fullfile(woi1_dir,filename));
    pial_woi1_trials(:,end+1)=trial_mesh.cdata(:);
end

% Load all pial data from woi2
pial_woi2_trials=[];
[files,dirs] = spm_select('List', woi2_dir, ['^pial_r' subj_info.subj_id '_' num2str(session_num) '_1_t' num2str(woi2(1)) '_' num2str(woi2(2)) '_f' num2str(foi(1)) '_' num2str(foi(2)) '.*\.gii']);
for f=1:size(files,1)
    filename=deblank(files(f,:));
    trial_mesh=gifti(fullfile(woi2_dir,filename));
    pial_woi2_trials(:,end+1)=trial_mesh.cdata(:);
end
pial_diff=pial_woi1_trials-pial_woi2_trials;

% Load all white matter data from woi1
white_woi1_trials=[];
[files,dirs] = spm_select('List', woi1_dir, ['^white_r' subj_info.subj_id '_' num2str(session_num) '_1_t' num2str(woi1(1)) '_' num2str(woi1(2)) '_f' num2str(foi(1)) '_' num2str(foi(2)) '.*\.gii']);
for f=1:size(files,1)
    filename=deblank(files(f,:));
    trial_mesh=gifti(fullfile(woi1_dir,filename));
    white_woi1_trials(:,end+1)=trial_mesh.cdata(:);
end

% Load all white matter data from woi2
white_woi2_trials=[];
[files,dirs] = spm_select('List', woi2_dir, ['^white_r' subj_info.subj_id '_' num2str(session_num) '_1_t' num2str(woi2(1)) '_' num2str(woi2(2)) '_f' num2str(foi(1)) '_' num2str(foi(2)) '.*\.gii']);
for f=1:size(files,1)
    filename=deblank(files(f,:));
    trial_mesh=gifti(fullfile(woi2_dir, filename));
    white_woi2_trials(:,end+1)=trial_mesh.cdata(:);
end
white_diff=white_woi1_trials-white_woi2_trials;

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
