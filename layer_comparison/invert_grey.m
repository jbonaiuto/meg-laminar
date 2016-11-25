function invert_grey(subj_info, session_num, zero_event, foi, woi, baseline, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'patch_size',0.4, 'surf_dir', '', 'mri_dir', '', 'init', true, 'coreg', true, 'invert', true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end
if length(params.surf_dir)==0
    params.surf_dir=fullfile(params.data_dir,'surf');
end
if length(params.mri_dir)==0
    params.mri_dir=fullfile(params.data_dir,'mri');
end

data_dir=fullfile(params.data_dir,'analysis', subj_info.subj_id, num2str(session_num));
data_file_name=fullfile(data_dir, sprintf('rc%s_Tafdf%d.mat', zero_event, session_num));

% Create directory for inversion results
foi_dir=fullfile(data_dir, 'grey_coreg', ['p' num2str(params.patch_size)], ['f' num2str(foi(1)) '_' num2str(foi(2))]);
if exist(foi_dir,'dir')~=7
    mkdir(foi_dir);
end
coreg_file_name=fullfile(foi_dir, sprintf('%s_%d.mat', subj_info.subj_id, session_num));
removed_file_name=fullfile(foi_dir, sprintf('r%s_%d.mat', subj_info.subj_id, session_num));
bc_file_name=fullfile(foi_dir, sprintf('br%s_%d.mat', subj_info.subj_id, session_num));

spm('defaults', 'EEG');
spm_jobman('initcfg'); 

if params.init
    % Copy file to foi_dir
    clear jobs
    matlabbatch={};
    batch_idx=1;
    matlabbatch{batch_idx}.spm.meeg.other.copy.D = {data_file_name};
    matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = coreg_file_name;
    batch_idx=batch_idx+1;

    % Remove bad trials
    matlabbatch{batch_idx}.spm.meeg.preproc.remove.D = {coreg_file_name};
    matlabbatch{batch_idx}.spm.meeg.preproc.remove.prefix = 'r';
    batch_idx=batch_idx+1;

    % Set EEG channels to other
    matlabbatch{batch_idx}.spm.meeg.preproc.prepare.D = {removed_file_name};
    matlabbatch{batch_idx}.spm.meeg.preproc.prepare.task{1}.settype.channels{1}.type = 'EEG';
    matlabbatch{batch_idx}.spm.meeg.preproc.prepare.task{1}.settype.newtype = 'Other';
    batch_idx=batch_idx+1;

    %%%%%% BASELINE CORRECT %%%%%%%%
    matlabbatch{batch_idx}.spm.meeg.preproc.bc.D = {removed_file_name};
    matlabbatch{batch_idx}.spm.meeg.preproc.bc.timewin = baseline;
    matlabbatch{batch_idx}.spm.meeg.preproc.bc.prefix = 'b';
    batch_idx=batch_idx+1;
    
    spm_jobman('run', matlabbatch);

    % Relabel trials to all be same condition
    load(bc_file_name);
    D.condlist={zero_event};
    for trial_idx=1:length(D.trials)
        D.trials(trial_idx).label=zero_event;
    end
    save(bc_file_name,'D');
    
end

if params.coreg
    spm_jobman('initcfg'); 
    clear jobs
    matlabbatch={};
    batch_idx=1;

    % Coregister with surface
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.D = {bc_file_name};
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.comment = 'grey';
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {fullfile(params.mri_dir,[subj_info.subj_id subj_info.birth_date], [subj_info.headcast_t1 ',1'])};
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed-ds_pial.hires.deformed.surf.gii')};
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshes.custom.iskull = {''};
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshes.custom.oskull = {''};
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshes.custom.scalp = {''};
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshres = 2;
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = subj_info.nas;
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = subj_info.lpa;
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = subj_info.rpa;
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
    matlabbatch{batch_idx}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    spm_jobman('run',matlabbatch);
end

if params.invert
    % Run the inversion
    spm_jobman('initcfg'); 
    clear jobs
    matlabbatch={};
    batch_idx=1;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.D = {bc_file_name};
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.val = 1;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.whatconditions.all = 1;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.invtype = 'EBB';
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.woi = woi;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.foi = foi;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.hanning = 0;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.npatches = 512;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.niter = 1;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm = params.patch_size;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.mselect = 0;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.nsmodes = 180;
    %matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.umodes = 'C:\RAW_TDCS_data\Analysis_020715\GB\Ugeneric.mat';
    %matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = 0;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = 16;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.priors.priorsmask = {''};
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 1;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
    matlabbatch{batch_idx}.spm.meeg.source.invertiter.modality = {'MEG'};
    spm_jobman('run',matlabbatch);
end