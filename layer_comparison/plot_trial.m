function plot_trial(subj_info, data_dir, subj_surf_fname, baseline)

a=load(fullfile(data_dir,'times.mat'));

ntimes=length(a.times);
baseline_idx=intersect(find(a.times>=baseline(1)),find(a.times<=baseline(2)));
subj_surf=gifti(subj_surf_fname);
min_val=100000;
max_val=-100000;
baseline=zeros(size(subj_surf.vertices,1),length(baseline_idx));
for j=1:length(baseline_idx)
    pial_surf = gifti(fullfile(data_dir,sprintf('pialwhite_t%d.gii',baseline_idx(j))));
    baseline(:,j)=pial_surf.cdata(:);
end
mean_baseline=mean(baseline,2);
for j=1:ntimes
    pial_surf = gifti(fullfile(data_dir,sprintf('pialwhite_t%d.gii',j)));
    bc=(pial_surf.cdata(:)-mean_baseline)./mean_baseline;
    max_time_val=max(bc);
    if max_time_val>max_val
        max_val=max_time_val;
    end
    min_time_val=min(bc);
    if min_time_val<min_val
        min_val=min_time_val;
    end
end
for j=1:ntimes
    pial_surf = gifti(fullfile(data_dir,sprintf('pialwhite_t%d.gii',j)));
    pial_surf.cdata=(pial_surf.cdata(:)-mean_baseline)./mean_baseline;
    plot_surface_metric(subj_info, subj_surf, pial_surf, 'topdown', 'output_file', fullfile(data_dir, sprintf('pialwhite_t%d.png',j)), 'limits', [min_val max_val], 'title', sprintf('%.1fms',a.times(j)));
    if mod(j,10)==0
        close all;
    end
end

