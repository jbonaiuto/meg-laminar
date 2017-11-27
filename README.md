Analysis code for Bonaiuto, et al "Laminar-specific cortical dynamics in human visual and sensorimotor cortices"

Requirements:

* RMAOV2: https://uk.mathworks.com/matlabcentral/fileexchange/5578-rmaov2
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
