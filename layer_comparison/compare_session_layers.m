function compare_session_layers(subj_info, session_num, zero_event, foi, wois, baseline, comparison_names, varargin)

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

foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id, num2str(session_num), 'grey_coreg', params.inv_type, ['p' num2str(params.patch_size)], zero_event, ['f' num2str(foi(1)) '_' num2str(foi(2))]);

D=spm_eeg_load(fullfile(foi_dir,sprintf('r%s_%d.mat', subj_info.subj_id, session_num)));
goodchans=D.indchantype('MEGGRAD','good');
Dgood=squeeze(D(goodchans,:,:));
M=D.inv{1}.inverse.M;
U=D.inv{1}.inverse.U{1};
n_combined_vertices=size(M,1);
n_vertices=round(n_combined_vertices/2);
times=D.inv{1}.inverse.pst;
baseline_idx=intersect(find(times>=baseline(1)),find(times<=baseline(2)));
ntrials=size(Dgood,3);
ncomparisons=length(comparison_names);
woi_vals=zeros(ncomparisons,n_combined_vertices,ntrials);
baselines=zeros(n_combined_vertices,ntrials);
parfor t=1:ntrials
%for t=1:ntrials
    d1=squeeze(Dgood(:,:,t));
    Dtrial=M*U*d1;
    baselines(:,t)=sum(Dtrial(:,baseline_idx).^2,2);
    for i=1:ncomparisons
        woi_idx=intersect(find(times>=wois(i,1)),find(times<=wois(i,2)));        
        woi_vals(i,:,t)=sum(Dtrial(:,woi_idx).^2,2);
    end
end
mean_baseline=mean(baselines,2);

pial_diffs=zeros(ncomparisons,n_vertices,ntrials);
white_diffs=zeros(ncomparisons,n_vertices,ntrials);
for t=1:ntrials
    for i=1:ncomparisons
        diff=(squeeze(woi_vals(i,:,t))-mean_baseline')./mean_baseline';
        pial_diffs(i,:,t)=diff(n_vertices+1:end);
        white_diffs(i,:,t)=diff(1:n_vertices);
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
