% List of open inputs
% Inversion results: M/EEG datasets - cfg_files
% Inversion results: M/EEG datasets - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'D:\pred_coding\src\matlab\analysis\layer_comparison\inversion_results_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Inversion results: M/EEG datasets - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Inversion results: M/EEG datasets - cfg_files
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
