function invert_grey_coregerr(subj_info, session_num, contrast, idx, varargin)

% Parse inputs
defaults = struct('data_dir', 'd:/pred_coding/derivatives/spm12', 'inv_type', 'EBB',...
    'patch_size',0.4, 'surf_dir', 'd:/pred_coding/derivatives/freesurfer', 'mri_dir', 'd:/pred_coding/', 'init', true,...
    'coreg', true, 'invert', true, 'shift_magnitude', 10);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

data_dir=fullfile(params.data_dir,subj_info.subj_id, sprintf('ses-0%d',session_num));

% Create directory for inversion results
origfoi_dir=fullfile(data_dir, 'grey_coreg', params.inv_type, ['p' num2str(params.patch_size)],...
    contrast.zero_event, ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);
bc_file_name=fullfile(origfoi_dir, sprintf('br%s_%d.mat', subj_info.subj_id, session_num));

foi_dir=fullfile(data_dir, 'grey_coreg', params.inv_type, ['p' num2str(params.patch_size)],...
    contrast.zero_event, ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))],...
    'coregerr',num2str(params.shift_magnitude), num2str(idx));
if exist(foi_dir,'dir')~=7
    mkdir(foi_dir);
end
coregerr_bc_file_name=fullfile(foi_dir, sprintf('br%s_%d.mat', subj_info.subj_id, session_num));

spm('defaults', 'EEG');
spm_jobman('initcfg'); 


% Copy file to foi_dir
clear jobs
matlabbatch={};
batch_idx=1;
matlabbatch{batch_idx}.spm.meeg.other.copy.D = {bc_file_name};
matlabbatch{batch_idx}.spm.meeg.other.copy.outfile = coregerr_bc_file_name;
batch_idx=batch_idx+1;

% Coregister with surface
matlabbatch{batch_idx}.spm.meeg.source.headmodel.D = {coregerr_bc_file_name};
matlabbatch{batch_idx}.spm.meeg.source.headmodel.val = 1;
matlabbatch{batch_idx}.spm.meeg.source.headmodel.comment = 'grey';
matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {fullfile(params.mri_dir,subj_info.subj_id, 'anat', [subj_info.headcast_t1 ',1'])};
matlabbatch{batch_idx}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {fullfile(params.surf_dir,subj_info.subj_id,'ds_white.hires.deformed-ds_pial.hires.deformed.surf.gii')};
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
batch_idx=batch_idx+1;

% Run the inversion
matlabbatch{batch_idx}.spm.meeg.source.invertiter.D = {coregerr_bc_file_name};
matlabbatch{batch_idx}.spm.meeg.source.invertiter.val = 1;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.whatconditions.all = 1;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.invtype = params.inv_type;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.woi = contrast.invwoi;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.foi = contrast.foi;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.hanning = 0;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.npatches = 512;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.niter = 1;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm = params.patch_size;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.mselect = 0;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.nsmodes = 180;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = 16;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.priors.priorsmask = {''};
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 1;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.shuffle = 0;
matlabbatch{batch_idx}.spm.meeg.source.invertiter.modality = {'MEG'};
spm_jobman('run',matlabbatch);
