%-----------------------------------------------------------------------
% Job saved on 27-Oct-2016 12:13:59 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6470)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.stats.factorial_design.dir = '<UNDEFINED>';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'congruence';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'coherence';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 3;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1
                                                                    1];
%%
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans = '<UNDEFINED>';
%%
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [1
                                                                    2];
%%
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = '<UNDEFINED>';
%%
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [1
                                                                    3];
%%
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = '<UNDEFINED>';
%%
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [2
                                                                    1];
%%
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans = '<UNDEFINED>';
%%
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(5).levels = [2
                                                                    2];
%%
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(5).scans = '<UNDEFINED>';
%%
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(6).levels = [2
                                                                    3];
%%
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(6).scans = '<UNDEFINED>';
%%
matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Negative effect of condition';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [-0.167 -0.167 -0.167 -0.167 -0.167 -0.167];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = '<UNDEFINED>';
matlabbatch{4}.spm.stats.results.print = 'png';
matlabbatch{4}.spm.stats.results.write.none = 1;
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 5;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{5}.spm.stats.results.conspec.extent = 0;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = '<UNDEFINED>';
matlabbatch{5}.spm.stats.results.print = 'png';
matlabbatch{5}.spm.stats.results.write.none = 1;
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 11;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = '<UNDEFINED>';
matlabbatch{6}.spm.stats.results.print = 'png';
matlabbatch{6}.spm.stats.results.write.none = 1;
