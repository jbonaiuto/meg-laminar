function subject_sensor_tf(subj_info, zero_event, varargin)

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

matlabbatch{batch_idx}.spm.stats.factorial_design.dir = {fullfile(params.data_dir, 'analysis', ...
    subj_info.subj_id, [type '_rtf_rc' zero_event '_Tafdf'])};

for session_num=1:length(subj_info.sessions)
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.anovaw.fsubject(session_num).scans = {};

    session_conditions=[];
    % iterate through each run    
    for run_num=1:subj_info.sessions(session_num)        
        for cond_idx=1:length(conditions)
            img_path=fullfile(params.data_dir, 'analysis', ...
                subj_info.subj_id, num2str(session_num), ...
                [type '_rtf_rc' zero_event '_Tafdf' num2str(session_num)],...
                ['scondition_' conditions{cond_idx} '.nii']);        
            V=spm_vol(img_path);
            for j=1:length(V)
                matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(session_num).scans{end+1,1}=fullfile(sprintf('%s,%d', img_path, j));
                session_conditions(end+1)=cond_idx;
            end
        end
    end
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.anovaw.fsubject(session_num).conds = session_conditions';
end

matlabbatch{batch_idx}.spm.stats.factorial_design.des.anovaw.dept = 0;
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

matlabbatch{batch_idx}.spm.stats.fmri_est.spmmat(1) = {fullfile(params.data_dir, 'analysis', ...
    subj_info.subj_id, [type '_rtf_rc' zero_event '_Tafdf'],'SPM.mat')};
matlabbatch{batch_idx}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{batch_idx}.spm.stats.fmri_est.method.Classical = 1;
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = {fullfile(params.data_dir, 'analysis', ...
    subj_info.subj_id, [type '_rtf_rc' zero_event '_Tafdf'],'SPM.mat')};
matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.name = 'congruent - low';
matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.weights = [1 -0.2 -0.2 -0.2 -0.2 -0.2 zeros(1,length(subj_info.sessions))];
matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.consess{2}.tcon.name = 'congruent - med';
matlabbatch{batch_idx}.spm.stats.con.consess{2}.tcon.weights = [-0.2 1 -0.2 -0.2 -0.2 -0.2 zeros(1,length(subj_info.sessions))];
matlabbatch{batch_idx}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.consess{3}.tcon.name = 'congruent - high';
matlabbatch{batch_idx}.spm.stats.con.consess{3}.tcon.weights = [-0.2 -0.2 1 -0.2 -0.2 -0.2 zeros(1,length(subj_info.sessions))];
matlabbatch{batch_idx}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.consess{4}.tcon.name = 'incongruent - low';
matlabbatch{batch_idx}.spm.stats.con.consess{4}.tcon.weights = [-0.2 -0.2 -0.2 1 -0.2 -0.2 zeros(1,length(subj_info.sessions))];
matlabbatch{batch_idx}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.consess{5}.tcon.name = 'incongruent - med';
matlabbatch{batch_idx}.spm.stats.con.consess{5}.tcon.weights = [-0.2 -0.2 -0.2 -0.2 1 -0.2 zeros(1,length(subj_info.sessions))];
matlabbatch{batch_idx}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.consess{6}.tcon.name = 'incongruent - high';
matlabbatch{batch_idx}.spm.stats.con.consess{6}.tcon.weights = [-0.2 -0.2 -0.2 -0.2 -0.2 1 zeros(1,length(subj_info.sessions))];
matlabbatch{batch_idx}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.consess{7}.tcon.name = 'congruent > incongruent';
matlabbatch{batch_idx}.spm.stats.con.consess{7}.tcon.weights = [1 1 1 -1 -1 -1 zeros(1,length(subj_info.sessions))];
matlabbatch{batch_idx}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.consess{8}.tcon.name = 'congruent < incongruent';
matlabbatch{batch_idx}.spm.stats.con.consess{8}.tcon.weights = [-1 -1 -1 1 1 1 zeros(1,length(subj_info.sessions))];
matlabbatch{batch_idx}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.consess{9}.tcon.name = 'coherence - increasing';
matlabbatch{batch_idx}.spm.stats.con.consess{9}.tcon.weights = [-1 0 1 -1 0 1 zeros(1,length(subj_info.sessions))];
matlabbatch{batch_idx}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.consess{10}.tcon.name = 'coherence - decreasing';
matlabbatch{batch_idx}.spm.stats.con.consess{10}.tcon.weights = [1 0 -1 1 0 -1 zeros(1,length(subj_info.sessions))];
matlabbatch{batch_idx}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.consess{11}.tcon.name = 'prediction error';
matlabbatch{batch_idx}.spm.stats.con.consess{11}.tcon.weights = [1 0 -1 -1 0 1 zeros(1,length(subj_info.sessions))];
matlabbatch{batch_idx}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.consess{12}.tcon.name = 'inverse prediction error';
matlabbatch{batch_idx}.spm.stats.con.consess{12}.tcon.weights = [-1 0 1 1 0 -1 zeros(1,length(subj_info.sessions))];
matlabbatch{batch_idx}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.delete = 0;
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