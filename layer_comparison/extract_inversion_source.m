function extract_inversion_source(subj_info, session_num, contrast, foi_dir, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'inv_type', 'EBB',...
    'patch_size',0.4,'surf_dir','d:/pred_coding/surf');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

coreg_file_name=fullfile(foi_dir, sprintf('br%s_%d.mat', subj_info.subj_id, session_num));

% Extract source at woi1 and woi2
spm_jobman('initcfg'); 
clear jobs
matlabbatch={};
batch_idx=1;

wois=[contrast.comparison_woi; contrast.baseline_woi];
for w=1:size(wois,1)
    woi=wois(w,:);
    matlabbatch{batch_idx}.spm.meeg.source.results.D = {coreg_file_name};
    matlabbatch{batch_idx}.spm.meeg.source.results.val = 1;
    matlabbatch{batch_idx}.spm.meeg.source.results.woi = woi;
    matlabbatch{batch_idx}.spm.meeg.source.results.foi = contrast.foi;
    matlabbatch{batch_idx}.spm.meeg.source.results.ctype = 'trials';
    matlabbatch{batch_idx}.spm.meeg.source.results.space = 0;
    matlabbatch{batch_idx}.spm.meeg.source.results.format = 'mesh';
    matlabbatch{batch_idx}.spm.meeg.source.results.smoothing = 8;
    batch_idx=batch_idx+1;
end
spm_jobman('run',matlabbatch);

for w=1:size(wois,1)
    woi=wois(w,:);
    
    % Save pial comparison
    woi_dir=fullfile(foi_dir, ['t' num2str(woi(1)) '_' num2str(woi(2))]);
    if exist(woi_dir,'dir')~=7
        mkdir(woi_dir);
    end
    delete(fullfile(woi_dir,'*'));
    movefile(fullfile(foi_dir, sprintf('br%s_%d_1_t%d_%d_f%d_%d_*', subj_info.subj_id, session_num, woi(1), woi(2), contrast.foi(1), contrast.foi(2))), woi_dir);
    
    % Split pial and grey sources
    split_inversion_results(woi_dir);
 
 end




