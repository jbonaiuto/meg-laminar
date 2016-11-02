%-----------------------------------------------------------------------
% Job saved on 02-Nov-2016 13:47:43 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6470)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.meeg.source.headmodel.D = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = 'grey';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.iskull = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.oskull = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.scalp = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
matlabbatch{2}.spm.meeg.source.invertiter.D(1) = cfg_dep('Head model specification: M/EEG dataset(s) with a forward model', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
matlabbatch{2}.spm.meeg.source.invertiter.val = 1;
matlabbatch{2}.spm.meeg.source.invertiter.whatconditions.all = 1;
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.invtype = 'EBB';
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.woi = '<UNDEFINED>';
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.foi = '<UNDEFINED>';
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.hanning = 0;
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.npatches = 512;
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.niter = 1;
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm = '<UNDEFINED>';
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.mselect = 0;
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.nsmodes = 180;
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.umodes = '';
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = 16;
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.priors.priorsmask = {''};
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 1;
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
matlabbatch{2}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
matlabbatch{2}.spm.meeg.source.invertiter.modality = {'MEG'};
