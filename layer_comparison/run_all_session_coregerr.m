function run_all_session_coregerr(subj_info, session_num, contrasts, varargin)

defaults = struct('data_dir', 'd:/pred_coding/derivatives/spm12', 'inv_type', 'EBB',...
    'patch_size',0.4, 'surf_dir', 'd:/pred_coding/derivatives/freesurfer', 'mri_dir', 'd:/pred_coding',... 
    'invert',true, 'extract', true, 'compare', true, 'iterations',10,...
    'shift_magnitude', 10 ,'shift_var', 2.5,...
    'rotation_magnitude', 10, 'rotation_var', 2.5);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

orig_nas=subj_info.nas;
orig_lpa=subj_info.lpa;
orig_rpa=subj_info.rpa;

rot_mag_rad=params.rotation_magnitude*pi/180.0;
rot_var_rad=params.rotation_var*pi/180.0;

for idx=1:params.iterations
    shift_vec=randn(1,3);
    shift_vec=shift_vec./sqrt(sum(shift_vec.^2)).*(randn*params.shift_var+params.shift_magnitude);
    rot_vec=randn(1,3);
    rot_vec=rot_vec./sqrt(sum(rot_vec.^2)).*(randn*rot_var_rad+rot_mag_rad);
    P=[shift_vec rot_vec];
    [A]=spm_matrix(P);
    origfid=[orig_nas; orig_lpa; orig_rpa];
    meanfid=mean(origfid);
    zeromeanfid=[origfid-repmat(meanfid,3,1) ones(3,1)];
    newfid=(A*zeromeanfid')';
    newfid=newfid(:,1:3)+repmat(meanfid,3,1);
    
    subj_info.nas=newfid(1,:);
    subj_info.lpa=newfid(2,:);
    subj_info.rpa=newfid(3,:);
    
    for i=1:length(contrasts)
        contrast=contrasts(i);
        run_session_foi_layer_comparison_coregerr(subj_info, session_num,...
            contrast, idx, 'data_dir', params.data_dir, 'inv_type', params.inv_type,...
            'patch_size',params.patch_size, 'surf_dir', params.surf_dir,...
            'mri_dir', params.mri_dir, 'invert',params.invert,'extract', params.extract,...
            'compare', params.compare, 'shift_magnitude', params.shift_magnitude);
    end        
end
