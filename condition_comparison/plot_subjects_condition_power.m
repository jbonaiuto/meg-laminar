function cond_results=plot_subjects_condition_power(subjects, contrast, varargin)

% Parse inputs
defaults = struct('data_dir','d:/pred_coding/derivatives/spm12',...
    'surf_dir', 'D:/pred_coding/derivatives/freesurfer',...
    'inv_type','EBB','patch_size',0.4,...
    'thresh_percentile',80,'roi_type','mean','recompute', false, 'recompute_roi',false,...
    'correct_only', true, 'plot_subjects', false, 'plot', true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

spm('defaults','eeg');

% Split trials by condition
conditions={'congruent - low','congruent - med','congruent - high',...
    'incongruent - low','incongruent - med','incongruent - high'};
coherence_conditions={'low','med','high'};
congruence_conditions={'congruent','incongruent'};
accuracy_conditions={'correct', 'incorrect'};


cond_results.pial_trials_woi=dict();
cond_results.wm_trials_woi=dict();
cond_results.norm_pial_trials_woi=dict();
cond_results.norm_wm_trials_woi=dict();

subj_ids={};

for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    subj_ids{subj_idx}=subj_info.subj_id;    
    
    subj_cond_results=plot_subject_condition_power(subj_info,...
        contrast, 'data_dir', params.data_dir,'surf_dir', params.surf_dir,...
        'inv_type', params.inv_type, 'patch_size', params.patch_size,...
        'thresh_percentile', params.thresh_percentile,...
        'roi_type', params.roi_type, 'recompute', params.recompute, 'recompute_roi',params.recompute_roi,...
        'correct_only', params.correct_only, 'plot', params.plot_subjects);
    
    pial_scale_factor=get_scale_factor_woi(accuracy_conditions, subj_cond_results.pial_trials_woi, contrast);
    wm_scale_factor=get_scale_factor_woi(accuracy_conditions, subj_cond_results.wm_trials_woi, contrast);
    for cond_idx=1:length(accuracy_conditions)
        condition=accuracy_conditions{cond_idx};
        pial_woi=mean(subj_cond_results.pial_trials_woi(condition));
        wm_woi=mean(subj_cond_results.wm_trials_woi(condition));
        
        norm_pial_woi=pial_woi./pial_scale_factor;
        norm_wm_woi=wm_woi./wm_scale_factor;
        cond_results.pial_trials_woi(condition)=[cond_results.pial_trials_woi(condition) pial_woi];
        cond_results.wm_trials_woi(condition)=[cond_results.wm_trials_woi(condition) wm_woi];
        cond_results.norm_pial_trials_woi(condition)=[cond_results.norm_pial_trials_woi(condition) norm_pial_woi];
        cond_results.norm_wm_trials_woi(condition)=[cond_results.norm_wm_trials_woi(condition) norm_wm_woi];
    end
    
    x_woi=dict();
    x_woi('all-pial')=subj_cond_results.pial_trials_woi('all');
    x_woi('all-white')=subj_cond_results.wm_trials_woi('all');
    scale_factor=get_scale_factor_woi({'all-pial','all-white'}, x_woi, contrast);
    cond_pial_woi=mean(subj_cond_results.pial_trials_woi('all'));
    cond_wm_woi=mean(subj_cond_results.wm_trials_woi('all'));
    
    norm_cond_pial_woi=cond_pial_woi./scale_factor;
    norm_cond_wm_woi=cond_wm_woi./scale_factor;
    cond_results.pial_trials_woi('all')=[cond_results.pial_trials_woi('all') cond_pial_woi];
    cond_results.wm_trials_woi('all')=[cond_results.wm_trials_woi('all') cond_wm_woi];
    cond_results.norm_pial_trials_woi('all')=[cond_results.norm_pial_trials_woi('all') norm_cond_pial_woi];
    cond_results.norm_wm_trials_woi('all')=[cond_results.norm_wm_trials_woi('all') norm_cond_wm_woi];
    
    pial_scale_factor=get_scale_factor_woi(conditions, subj_cond_results.pial_trials_woi, contrast);
    wm_scale_factor=get_scale_factor_woi(conditions, subj_cond_results.wm_trials_woi, contrast);
    for cond_idx=1:length(conditions)
        condition=conditions{cond_idx};
        cond_pial_woi=mean(subj_cond_results.pial_trials_woi(condition));
        cond_wm_woi=mean(subj_cond_results.wm_trials_woi(condition));
        if strcmp(contrast.direction,'positive') && abs(cond_pial_woi)>abs(pial_scale_factor) && cond_pial_woi<pial_scale_factor
            norm_cond_pial_woi=(pial_scale_factor-(cond_pial_woi-pial_scale_factor))./pial_scale_factor;
        else
            norm_cond_pial_woi=cond_pial_woi./pial_scale_factor;
        end
        if strcmp(contrast.direction,'positive') && abs(cond_wm_woi)>abs(wm_scale_factor) && cond_wm_woi<wm_scale_factor
            norm_cond_wm_woi=(wm_scale_factor-(cond_wm_woi-wm_scale_factor))./wm_scale_factor;
        else
            norm_cond_wm_woi=cond_wm_woi./wm_scale_factor;
        end
        
        cond_results.pial_trials_woi(condition)=[cond_results.pial_trials_woi(condition) cond_pial_woi];
        cond_results.wm_trials_woi(condition)=[cond_results.wm_trials_woi(condition) cond_wm_woi];
        cond_results.norm_pial_trials_woi(condition)=[cond_results.norm_pial_trials_woi(condition) norm_cond_pial_woi];
        cond_results.norm_wm_trials_woi(condition)=[cond_results.norm_wm_trials_woi(condition) norm_cond_wm_woi];
    end
    
    pial_scale_factor=get_scale_factor_woi(coherence_conditions, subj_cond_results.pial_trials_woi, contrast);
    wm_scale_factor=get_scale_factor_woi(coherence_conditions, subj_cond_results.wm_trials_woi, contrast);
    for cond_idx=1:length(coherence_conditions)
        % Get trials in this condition - congruent and incongruent
        condition=coherence_conditions{cond_idx};
        cond_pial_woi=mean(subj_cond_results.pial_trials_woi(condition));
        cond_wm_woi=mean(subj_cond_results.wm_trials_woi(condition));
        if strcmp(contrast.direction,'positive') && abs(cond_pial_woi)>abs(pial_scale_factor) && cond_pial_woi<pial_scale_factor
            norm_cond_pial_woi=(pial_scale_factor-(cond_pial_woi-pial_scale_factor))./pial_scale_factor;
        else
            norm_cond_pial_woi=cond_pial_woi./pial_scale_factor;
        end
        if strcmp(contrast.direction,'positive') && abs(cond_wm_woi)>abs(wm_scale_factor) && cond_wm_woi<wm_scale_factor
            norm_cond_wm_woi=(wm_scale_factor-(cond_wm_woi-wm_scale_factor))./wm_scale_factor;
        else
            norm_cond_wm_woi=cond_wm_woi./wm_scale_factor;
        end
        
        cond_results.pial_trials_woi(condition)=[cond_results.pial_trials_woi(condition) cond_pial_woi];
        cond_results.wm_trials_woi(condition)=[cond_results.wm_trials_woi(condition) cond_wm_woi];
        cond_results.norm_pial_trials_woi(condition)=[cond_results.norm_pial_trials_woi(condition) norm_cond_pial_woi];
        cond_results.norm_wm_trials_woi(condition)=[cond_results.norm_wm_trials_woi(condition) norm_cond_wm_woi];
    end
    
    pial_scale_factor=get_scale_factor_woi(congruence_conditions, subj_cond_results.pial_trials_woi, contrast);
    wm_scale_factor=get_scale_factor_woi(congruence_conditions, subj_cond_results.wm_trials_woi, contrast);
    for cond_idx=1:length(congruence_conditions)
        condition=congruence_conditions{cond_idx};
        cond_pial_woi=mean(subj_cond_results.pial_trials_woi(condition));
        cond_wm_woi=mean(subj_cond_results.wm_trials_woi(condition));
        if strcmp(contrast.direction,'positive') && abs(cond_pial_woi)>abs(pial_scale_factor) && cond_pial_woi<pial_scale_factor
            norm_cond_pial_woi=(pial_scale_factor-(cond_pial_woi-pial_scale_factor))./pial_scale_factor;
        else
            norm_cond_pial_woi=cond_pial_woi./pial_scale_factor;
        end
        if strcmp(contrast.direction,'positive') && abs(cond_wm_woi)>abs(wm_scale_factor) && cond_wm_woi<wm_scale_factor
            norm_cond_wm_woi=(wm_scale_factor-(cond_wm_woi-wm_scale_factor))./wm_scale_factor;
        else
            norm_cond_wm_woi=cond_wm_woi./wm_scale_factor;
        end
        
        cond_results.pial_trials_woi(condition)=[cond_results.pial_trials_woi(condition) cond_pial_woi];
        cond_results.wm_trials_woi(condition)=[cond_results.wm_trials_woi(condition) cond_wm_woi];
        cond_results.norm_pial_trials_woi(condition)=[cond_results.norm_pial_trials_woi(condition) norm_cond_pial_woi];
        cond_results.norm_wm_trials_woi(condition)=[cond_results.norm_wm_trials_woi(condition) norm_cond_wm_woi];
    end
end

x_woi=dict();
x_woi('all-pial')=cond_results.norm_pial_trials_woi('all');
x_woi('all-white')=cond_results.norm_wm_trials_woi('all');

if params.plot
    
    fig=plot_power_woi(cond_results.norm_pial_trials_woi, accuracy_conditions, subj_ids);
    
    fig=plot_power_woi(cond_results.norm_wm_trials_woi, accuracy_conditions,subj_ids);
    
    fig=figure('position',[1 1 1185 950]);
    ax=subplot(2,2,1);
    plot_power_woi(x_woi, {'all-pial','all-white'}, subj_ids,'ax', ax);
    ax=subplot(2,2,2);
    plot_power_woi(cond_results.norm_pial_trials_woi, coherence_conditions, subj_ids,'ax', ax);
    title('pial');
    ax=subplot(2,2,3);
    plot_power_woi(cond_results.norm_pial_trials_woi, congruence_conditions,subj_ids,'ax', ax);
    title('pial');
    ax=subplot(2,2,4);
    plot_power_woi(cond_results.norm_pial_trials_woi, conditions, subj_ids,'ax', ax);
    title('pial');
    
    fig=figure('position',[1 1 1185 950]);
    ax=subplot(2,2,1);
    plot_power_woi(x_woi, {'all-pial','all-white'}, subj_ids,'ax', ax);
    ax=subplot(2,2,2);
    plot_power_woi(cond_results.norm_wm_trials_woi, coherence_conditions, subj_ids,'ax', ax);
    title('white');
    ax=subplot(2,2,3);
    plot_power_woi(cond_results.norm_wm_trials_woi, congruence_conditions, subj_ids,'ax', ax);
    title('white');
    ax=subplot(2,2,4);
    plot_power_woi(cond_results.norm_wm_trials_woi, conditions, subj_ids,'ax', ax);
    title('white');     
end
