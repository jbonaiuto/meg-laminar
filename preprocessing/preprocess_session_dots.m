function preprocess_session_dots(subj_info, session_num, varargin)

defaults = struct('data_dir', '/data/pred_coding', ...
    'url_prefix', 'http://fortressofjollitude.zapto.org/', 'dots_epoch',[-3500 0]);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Directory containing session data
session_dir=fullfile(params.data_dir,'scanning', subj_info.subj_id, num2str(session_num));
% Directory to put results
analysis_dir=fullfile(params.data_dir,'analysis',subj_info.subj_id,num2str(session_num));
if exist(analysis_dir,'dir')~=7
    mkdir(analysis_dir);
end

spm('defaults','eeg');
spm_jobman('initcfg');

clear jobs
matlabbatch={};
batch_idx=1;

matlabbatch{batch_idx}.spm.meeg.preproc.crop.D = {fullfile(analysis_dir, sprintf('rcinstr_Tafdf%d.mat',session_num))};
matlabbatch{batch_idx}.spm.meeg.preproc.crop.timewin = params.dots_epoch;
matlabbatch{batch_idx}.spm.meeg.preproc.crop.freqwin = [-Inf Inf];
matlabbatch{batch_idx}.spm.meeg.preproc.crop.channels{1}.all = 'all';
matlabbatch{batch_idx}.spm.meeg.preproc.crop.prefix = 'p';    
  
spm_jobman('run',matlabbatch);

clear jobs
matlabbatch={};
batch_idx=1;
    
matlabbatch{batch_idx}.spm.meeg.other.copy.D = {fullfile(analysis_dir, sprintf('prcinstr_Tafdf%d.mat',session_num))};
matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = fullfile(analysis_dir, sprintf('rcdots_Tafdf%d.mat', session_num));
batch_idx=batch_idx+1;

matlabbatch{batch_idx}.spm.meeg.other.delete.D = {fullfile(analysis_dir, sprintf('prcinstr_Tafdf%d.mat',session_num))};
batch_idx=batch_idx+1;

% Average
matlabbatch{batch_idx}.spm.meeg.averaging.average.D = {fullfile(analysis_dir, sprintf('rcdots_Tafdf%d.mat', session_num))};
matlabbatch{batch_idx}.spm.meeg.averaging.average.userobust.standard = false;
matlabbatch{batch_idx}.spm.meeg.averaging.average.plv = false;
matlabbatch{batch_idx}.spm.meeg.averaging.average.prefix = fullfile(analysis_dir, 'm');
batch_idx=batch_idx+1;

spm_jobman('run',matlabbatch);

