function group_stats(subjects, epoch_name, zero_evt)

spm('defaults', 'EEG');
spm_jobman('initcfg');
clear jobs;
matlabbatch={};
batch_idx=1;

dir=fullfile('C:\pred_coding\analysis\sensor_tf\stats',epoch_name);
mkdir(dir);

matlabbatch{batch_idx}.spm.stats.factorial_design.dir = {dir};

for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    
    scans={};
    sessions=[1:length(subj_info.sessions)];
    if strcmp(subj_info.subj_id,'nc')
        sessions=[3];
    end
    for session_idx=1:length(sessions)
        session_num=sessions(session_idx);
        fname=fullfile('C:\pred_coding\analysis',subj_info.subj_id,num2str(session_num),sprintf('r%s_tf_ffrc%s_Tafdf%d', epoch_name, zero_evt, session_num),'scondition_Undefined.nii');
        x=spm_vol(fname);
        for i=1:length(x)
            scans{end+1,1}=sprintf('%s,%d', fname, i);
        end
    end
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.anovaw.fsubject(subj_idx).scans = scans;
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.anovaw.fsubject(subj_idx).conds = ones(1,size(scans,1));
end
matlabbatch{batch_idx}.spm.stats.factorial_design.des.anovaw.dept = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.des.anovaw.variance = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.des.anovaw.gmsca = 0;
matlabbatch{batch_idx}.spm.stats.factorial_design.des.anovaw.ancova = 0;
matlabbatch{batch_idx}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{batch_idx}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{batch_idx}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{batch_idx}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.glonorm = 1;
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{batch_idx}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{batch_idx}.spm.stats.fmri_est.method.Classical = 1;
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.name = 'all';
matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.weights = [zeros(length(subjects),1) eye(length(subjects),length(subjects))-1/length(subjects)];
matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.delete = 1;
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{batch_idx}.spm.stats.results.spmmat = {fullfile(dir, 'SPM.mat')};
matlabbatch{batch_idx}.spm.stats.results.conspec.titlestr = '';
matlabbatch{batch_idx}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{batch_idx}.spm.stats.results.conspec.threshdesc = 'FDR';
matlabbatch{batch_idx}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{batch_idx}.spm.stats.results.conspec.extent = 500;
matlabbatch{batch_idx}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{batch_idx}.spm.stats.results.units = 4;
matlabbatch{batch_idx}.spm.stats.results.print = 'csv';
matlabbatch{batch_idx}.spm.stats.results.write.binary.basename = 'F_FDR';

spm_jobman('run', matlabbatch);

X=spm_vol(fullfile(dir, 'spmF_0001_F_FWE.nii'));
max_dim=max(X.dim);
coords=X.mat*[[1:max_dim]' [1:max_dim]' ones(max_dim,1) ones(max_dim,1)]';
freq_idx=[1:X.dim(1)];
times=coords(2,1:X.dim(2));
freqs=coords(1,freq_idx);
        
img=spm_read_vols(X);
img=squeeze(img(freq_idx,:,1,:));
alpha_data=ones(size(img));
alpha_data(find(img==1))=0;

fig=figure();
imagesc(times, freqs, img,'AlphaData',alpha_data);
set(gca,'ydir','normal');
xlabel('Time');
ylabel('Freq');
colorbar();
saveas(fig, fullfile(dir, sprintf('%s_mask.png',epoch_name)), 'png');