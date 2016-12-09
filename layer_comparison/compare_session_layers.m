function compare_session_layers(subj_info, session_num, zero_event, foi, wois, baseline_zero_event, baseline_woi, comparison_names, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'save_results', true, 'inv_type','EBB','patch_size', 0.4, 'white_pial_map', [],'surf_dir','');  %define default values
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

baseline_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id, num2str(session_num), 'grey_coreg', params.inv_type, ['p' num2str(params.patch_size)], baseline_zero_event, ['f' num2str(foi(1)) '_' num2str(foi(2))],['t' num2str(baseline_woi(1)) '_' num2str(baseline_woi(2))]);
[files,dirs] = spm_select('List', baseline_dir, '^pial.*\.gii');
ncomparisons=length(comparison_names);
ntrials=size(files,1);
nvertices=length(params.white_pial_map);
pial_baseline_vals=zeros(nvertices,ntrials);
white_baseline_vals=zeros(nvertices,ntrials);
for t=1:ntrials
    file_part=sprintf('_r%s_%d_1_t%d_%d_f%d_%d_1_%d',subj_info.subj_id,session_num,baseline_woi(1),baseline_woi(2),foi(1),foi(2),t);
    baseline_pial=gifti(fullfile(baseline_dir,sprintf('pial%s.gii',file_part)));
    pial_baseline_vals(:,t)=baseline_pial.cdata(:);
    baseline_white=gifti(fullfile(baseline_dir,sprintf('white%s.gii',file_part)));
    white_baseline_vals(:,t)=baseline_white.cdata(:);
end

foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id, num2str(session_num), 'grey_coreg', params.inv_type, ['p' num2str(params.patch_size)], zero_event, ['f' num2str(foi(1)) '_' num2str(foi(2))]);

pial_diffs=zeros(ncomparisons,nvertices,ntrials);
white_diffs=zeros(ncomparisons,nvertices,ntrials);
for i=1:ncomparisons
    woi_dir=fullfile(foi_dir,['t' num2str(wois(i,1)) '_' num2str(wois(i,2))]);
    for t=1:ntrials
        file_part=sprintf('_r%s_%d_1_t%d_%d_f%d_%d_1_%d',subj_info.subj_id,session_num,wois(i,1),wois(i,2),foi(1),foi(2),t);
        woi_pial=gifti(fullfile(woi_dir,sprintf('pial%s.gii',file_part)));
        %pial_diffs(i,:,t)=(woi_pial.cdata(:)-pial_baseline_vals(:,t))./pial_baseline_vals(:,t);
        pial_diffs(i,:,t)=woi_pial.cdata(:)-pial_baseline_vals(:,t);
        woi_white=gifti(fullfile(woi_dir,sprintf('white%s.gii',file_part)));
        %white_diffs(i,:,t)=(woi_white.cdata(:)-white_baseline_vals(:,t))./white_baseline_vals(:,t);
        white_diffs(i,:,t)=woi_white.cdata(:)-white_baseline_vals(:,t);
    end
end

if params.save_results
    for i=1:ncomparisons
        comparison_name=comparison_names{i};
        pial_diff=squeeze(pial_diffs(i,:,:));
        white_diff=squeeze(white_diffs(i,:,:));
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

end
