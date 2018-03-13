Analysis code for Bonaiuto, et al "Laminar-specific cortical dynamics in human visual and sensorimotor cortices"
https://www.biorxiv.org/content/early/2017/11/30/226274

Requirements:

* RMAOV2: https://uk.mathworks.com/matlabcentral/fileexchange/5578-rmaov2
* shadedErrorBar: https://fr.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar?
* barwitherr: https://fr.mathworks.com/matlabcentral/fileexchange/30639-barwitherr-errors-varargin-?focused=3845794&tab=function
* dict: http://uk.mathworks.com/matlabcentral/fileexchange/19381-lookuptable
* MEGSurfer: https://github.com/jbonaiuto/MEGsurfer

% To run analysis


% Create subject structure

    subjects=create_subject_structure();

% Create contrast structure

    contrasts=create_contrast_structure();


%%% Analyze subject behavior

    subjects_behavior(subjects);


%%% Reproducibility

% Extract fiducial coil coordinates

    extract_fiducials(subjects);

% Analyze fiducial coil movement

    analyze_fiducials(subjects);


% Sensor cluster TF

    run_tfs(subjects(3));
    extract_tfs(subjects(3));
    plot_tfs(subjects(3));


% ERPs

    plot_erps(subjects(3));


% Field maps

    plot_field_maps(subjects(3));


%%% Laminar analysis

% Run inversions

    for idx=1:length(subjects)
        for session_num=1:length(subjects(idx).sessions)
            run_all_session(subjects(idx), session_num, contrasts);
        end
        run_all_subject(subjects(idx), contrasts);
    end


% Plot example subjects

    for idx=1:length(contrasts)
        plot_classification_results_subject(subjects, contrasts(idx), 'thresh_percentile', 0, 'whole_brain', true, 'plot_ext', 'global_');
        plot_classification_results_subject(subjects, contrasts(idx), 'thresh_percentile', 80, 'whole_brain', true, 'plot_ext', 'func_');
        plot_classification_results_subject(subjects, contrasts(idx), 'thresh_percentile', 80, 'whole_brain', false, 'plot_ext', 'anat_');
    end


% Plot subject-level results

    for idx=1:length(contrasts)
        plot_classification_results_subjects(subjects, contrasts(idx));
    end


% Shuffled control

    for idx=1:length(subjects)
        for session_num=1:length(subjects(idx).sessions)
            run_all_session_shuffled(subjects(idx), session_num, contrasts);
        end
        run_all_subject_shuffled(subjects(idx), contrasts);
    end
    for idx=1:length(contrasts)
        plot_classification_results_subjects_shuffled(subjects, contrasts(idx));
    end



% Co-registration error control

    for idx=1:length(subjects)
        for session_num=1:length(subjects(idx).sessions)
            run_all_session_coregerr(subjects(idx), session_num, contrasts);
        end
        run_all_subject_coregerr(subjects(idx), contrasts);
    end
    for idx=1:length(contrasts)
        plot_classification_results_subjects_coregerr(subjects, contrasts(idx));
    end



% SNR control

    run_all_session_noise(subjects, contrasts);


%%% Condition comparison

    for idx=1:length(contrasts)
        plot_subjects_condition_power(subjects, contrasts(idx));
    end
