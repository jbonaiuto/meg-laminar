%-----------------------------------------------------------------------
% Job saved on 24-Oct-2016 14:18:24 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6470)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.meeg.tf.tf.D = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.tf.tf.channels{1}.type = 'MEG';
matlabbatch{1}.spm.meeg.tf.tf.frequencies = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100];
matlabbatch{1}.spm.meeg.tf.tf.timewin = [-Inf Inf];
matlabbatch{1}.spm.meeg.tf.tf.method.morlet.ncycles = 7;
matlabbatch{1}.spm.meeg.tf.tf.method.morlet.timeres = 0;
matlabbatch{1}.spm.meeg.tf.tf.method.morlet.subsample = 1;
matlabbatch{1}.spm.meeg.tf.tf.phase = 0;
matlabbatch{1}.spm.meeg.tf.tf.prefix = '';
matlabbatch{2}.spm.meeg.tf.rescale.D(1) = cfg_dep('Time-frequency analysis: M/EEG time-frequency power dataset', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dtfname'));
matlabbatch{2}.spm.meeg.tf.rescale.method.LogR.baseline.timewin = '<UNDEFINED>';
matlabbatch{2}.spm.meeg.tf.rescale.method.LogR.baseline.Db = [];
matlabbatch{2}.spm.meeg.tf.rescale.prefix = 'r';
matlabbatch{3}.spm.meeg.other.delete.D(1) = cfg_dep('Time-frequency analysis: M/EEG time-frequency power dataset', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dtfname'));
matlabbatch{4}.spm.meeg.images.convert2images.D(1) = cfg_dep('Time-frequency rescale: Rescaled TF Datafile', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{4}.spm.meeg.images.convert2images.mode = 'scalp x time';
matlabbatch{4}.spm.meeg.images.convert2images.conditions = cell(1, 0);
matlabbatch{4}.spm.meeg.images.convert2images.channels{1}.type = 'MEG';
matlabbatch{4}.spm.meeg.images.convert2images.timewin = '<UNDEFINED>';
matlabbatch{4}.spm.meeg.images.convert2images.freqwin = [2 100];
matlabbatch{4}.spm.meeg.images.convert2images.prefix = 'scalp_time_';
matlabbatch{5}.spm.meeg.images.convert2images.D(1) = cfg_dep('Time-frequency rescale: Rescaled TF Datafile', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{5}.spm.meeg.images.convert2images.mode = 'scalp x frequency';
matlabbatch{5}.spm.meeg.images.convert2images.conditions = cell(1, 0);
matlabbatch{5}.spm.meeg.images.convert2images.channels{1}.type = 'MEG';
matlabbatch{5}.spm.meeg.images.convert2images.timewin = '<UNDEFINED>';
matlabbatch{5}.spm.meeg.images.convert2images.freqwin = [2 100];
matlabbatch{5}.spm.meeg.images.convert2images.prefix = 'scalp_freq_';
matlabbatch{6}.spm.meeg.images.convert2images.D(1) = cfg_dep('Time-frequency rescale: Rescaled TF Datafile', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{6}.spm.meeg.images.convert2images.mode = 'time x frequency';
matlabbatch{6}.spm.meeg.images.convert2images.conditions = cell(1, 0);
matlabbatch{6}.spm.meeg.images.convert2images.channels{1}.type = 'MEG';
matlabbatch{6}.spm.meeg.images.convert2images.timewin = '<UNDEFINED>';
matlabbatch{6}.spm.meeg.images.convert2images.freqwin = [2 100];
matlabbatch{6}.spm.meeg.images.convert2images.prefix = 'time_freq_';
matlabbatch{7}.spm.meeg.other.delete.D(1) = cfg_dep('Time-frequency rescale: Rescaled TF Datafile', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
