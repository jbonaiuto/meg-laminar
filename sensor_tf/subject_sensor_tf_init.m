function subject_sensor_tf_init(subj_info, zero_event, varargin)

defaults = struct('data_dir', '/data/pred_coding');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults', 'EEG');
spm_jobman('initcfg');

conditions={'congruent-low','congruent-med','congruent-high',...
    'incongruent-low','incongruent-med','incongruent-high'};

run_analysis('scalp_freq',3);
run_analysis('scalp_time',2);
run_analysis('time_freq',4);


function run_analysis(type, units)
curr_dir=pwd;
spm_jobman('initcfg');
clear jobs;    
matlabbatch={};
batch_idx=1;

delete(fullfile(params.data_dir, 'analysis', ...
        subj_info.subj_id, [type '_rtf_rc' zero_event '_Tafdf'],'*'));

matlabbatch{batch_idx}.spm.stats.factorial_design.dir = {fullfile(params.data_dir, 'analysis', ...
        subj_info.subj_id, [type '_rtf_rc' zero_event '_Tafdf'])};
matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fac(1).name = 'subject';
matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fac(1).dept = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fac(1).variance = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fac(2).name = 'condition';
matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fac(2).dept = 0;
matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fac(2).variance = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;

for session_num=1:length(subj_info.sessions)
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(session_num).scans = {};
    
    % iterate through each run    
    for run_num=1:subj_info.sessions(session_num)        
        for cond_idx=1:length(conditions)
            img_path=fullfile(params.data_dir, 'analysis', ...
                subj_info.subj_id, num2str(session_num), ...
                [type '_rtf_rc' zero_event '_Tafdf' num2str(session_num)],...
                ['scondition_' conditions{cond_idx} '.nii']);        
            V=spm_vol(img_path);
            for j=1:length(V)
                matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(session_num).scans{end+1,1}=fullfile(sprintf('%s,%d', img_path, j));
            end
        end
    end
    nscans=length(matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(session_num).scans);
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(session_num).conds=ones(nscans,1);
end

matlabbatch{batch_idx}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 2;
matlabbatch{batch_idx}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{batch_idx}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{batch_idx}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{batch_idx}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.glonorm = 1;
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.stats.fmri_est.spmmat(1) = {fullfile(params.data_dir, 'analysis', ...
    subj_info.subj_id, [type '_rtf_rc' zero_event '_Tafdf'],'SPM.mat')};
matlabbatch{batch_idx}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{batch_idx}.spm.stats.fmri_est.method.Classical = 1;
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = {fullfile(params.data_dir, 'analysis', ...
    subj_info.subj_id, [type '_rtf_rc' zero_event '_Tafdf'],'SPM.mat')};
matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.name = 'positive';
matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.consess{2}.tcon.name = 'negative';
matlabbatch{batch_idx}.spm.stats.con.consess{2}.tcon.weights = -1;
matlabbatch{batch_idx}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.delete = 1;
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.stats.results.spmmat(1) = {fullfile(params.data_dir, 'analysis', ...
    subj_info.subj_id, [type '_rtf_rc' zero_event '_Tafdf'],'SPM.mat')};
matlabbatch{batch_idx}.spm.stats.results.conspec.titlestr = '';
matlabbatch{batch_idx}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{batch_idx}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{batch_idx}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{batch_idx}.spm.stats.results.conspec.extent = 0;
matlabbatch{batch_idx}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{batch_idx}.spm.stats.results.units = units;
matlabbatch{batch_idx}.spm.stats.results.print = 'png';
matlabbatch{batch_idx}.spm.stats.results.write.none = 1;
batch_idx=batch_idx+1;

spm_jobman('run', matlabbatch);

cd(curr_dir);
end
end