function plot_subjects_tfs(subjects, epoch_name, zero_evt, varargin)

defaults = struct('plot_subjects',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

mean_tfs=[];
for subj_idx=1:length(subjects)
    results=plot_subject_tfs(subjects(subj_idx),...
        epoch_name, zero_evt, 'plot',params.plot_subjects);
    mean_tfs(subj_idx,:,:)=results.session_mean_tfs;
end

fig=figure();
imagesc(results.times, results.freqs, squeeze(mean(mean_tfs)));
set(gca,'ydir','normal');
colorbar();
xlabel('Time');
ylabel('Freq');
