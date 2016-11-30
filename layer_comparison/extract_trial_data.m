function extract_trial_data(fname, outdir, varargin)

% Parse inputs
defaults = struct('max_trials',0,'parallel',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');
D=spm_eeg_load(fname);
goodchans=D.indchantype('MEGGRAD','good');
Dgood=squeeze(D(goodchans,:,:));
ntrials=size(Dgood,3);
if params.max_trials==0
    params.max_trials=ntrials;
end
M=D.inv{1}.inverse.M;
U=D.inv{1}.inverse.U{1};
times=D.inv{1}.inverse.pst;
ntimes=length(times);
nverts=size(M,1);
tc=zeros(nverts,ntimes);
if params.parallel
    parfor i=1:params.max_trials
        trial_tc=extract_individual_trial(i, M, U, Dgood, nverts, ntimes);
        tc=tc+trial_tc;
    end
else
    for i=1:params.max_trials
        disp(sprintf('trial %d',i));
        trial_tc=extract_individual_trial(i, M, U, Dgood, nverts, ntimes);
        tc=tc+trial_tc;
    end
end
tc=tc./params.max_trials;
for j=1:ntimes
    c=file_array(fullfile(outdir,sprintf('pialwhite_t%d.dat',j)),[nverts 1],'FLOAT32-LE',0,1,0);
    c(:)=tc(:,j);
    pial_surf = gifti;
    pial_surf.cdata = c;
    save(pial_surf, fullfile(outdir,sprintf('pialwhite_t%d.gii',j)), 'ExternalFileBinary');
end
save(fullfile(outdir,'times.mat'),'times');

