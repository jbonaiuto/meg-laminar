function analyze_fiducials(subjects)

spm('defaults','eeg');

subj_nas_session=[];
subj_lpa_session=[];
subj_rpa_session=[];

subj_labels={};
colors=[27,158,119;217,95,2;117,112,179;231,41,138;102,166,30;230,171,2;166,118,29;102,102,102]./255.0;

figure();
hold on;
for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    subj_labels{subj_idx}=subj_info.subj_id;

    session_init_nas=[];
    session_init_lpa=[];
    session_init_rpa=[];
    for session_num=1:length(subj_info.sessions)
        for run_num=1:subj_info.sessions(session_num)
            analysis_dir=fullfile('d:/pred_coding/derivatives/spm12',subj_info.subj_id, sprintf('ses-0%d',session_num));
            
            load(fullfile(analysis_dir, sprintf('mvmt_%d.mat',run_num)), 'mvmt');
                
            mvmt.nas_x=mvmt.nas_x(1:100:end);
            mvmt.nas_y=mvmt.nas_y(1:100:end);
            mvmt.nas_z=mvmt.nas_z(1:100:end);
            
            mvmt.lpa_x=mvmt.lpa_x(1:100:end);
            mvmt.lpa_y=mvmt.lpa_y(1:100:end);
            mvmt.lpa_z=mvmt.lpa_z(1:100:end);
            
            mvmt.rpa_x=mvmt.rpa_x(1:100:end);
            mvmt.rpa_y=mvmt.rpa_y(1:100:end);
            mvmt.rpa_z=mvmt.rpa_z(1:100:end);
            
            % If there is a jump of at least 0.1 in z direction, use up until jump
            end_idx=length(mvmt.nas_z);
            diff_z=abs(mvmt.nas_z(2:end)-mvmt.nas_z(1:end-1));
            if length(find(diff_z>0.001))
                end_idx=min(find(diff_z>0.001))-1;
            end

            run_nas_shift=[std(mvmt.nas_x(1:end_idx)) std(mvmt.nas_y(1:end_idx)) std(mvmt.nas_z(1:end_idx))];
            for dim_idx=1:3
                x=1+(dim_idx-2)*.25;
                plot(x,run_nas_shift(dim_idx).*1000,'o','MarkerEdgeColor',colors(subj_idx,:),'MarkerFaceColor',colors(subj_idx,:));
            end
            
            run_lpa_shift=[std(mvmt.lpa_x(1:end_idx)) std(mvmt.lpa_y(1:end_idx)) std(mvmt.lpa_z(1:end_idx))];
            for dim_idx=1:3
                x=2+(dim_idx-2)*.25;
                plot(x,run_lpa_shift(dim_idx).*1000,'o','MarkerEdgeColor',colors(subj_idx,:),'MarkerFaceColor',colors(subj_idx,:));
            end

            run_rpa_shift=[std(mvmt.rpa_x(1:end_idx)) std(mvmt.rpa_y(1:end_idx)) std(mvmt.rpa_z(1:end_idx))];
            for dim_idx=1:3
                x=3+(dim_idx-2)*.25;
                plot(x,run_rpa_shift(dim_idx).*1000,'o','MarkerEdgeColor',colors(subj_idx,:),'MarkerFaceColor',colors(subj_idx,:));
            end

            session_init_nas(end+1,:)=[mvmt.nas_x(1) mvmt.nas_y(1) mvmt.nas_z(1)];
            session_init_lpa(end+1,:)=[mvmt.lpa_x(1) mvmt.lpa_y(1) mvmt.lpa_z(1)];
            session_init_rpa(end+1,:)=[mvmt.rpa_x(1) mvmt.rpa_y(1) mvmt.rpa_z(1)];
        end
    end
    if size(session_init_nas,1)>1
        subj_nas_session(end+1,:)=std(session_init_nas);
        subj_lpa_session(end+1,:)=std(session_init_lpa);
        subj_rpa_session(end+1,:)=std(session_init_rpa);
    end
end
xlim([0.5 3.5]);
ylabel('SD (mm)');

figure();
hold on;
for dim_idx=1:size(subj_nas_session,2)
    x=1+(dim_idx-2)*.25;
    for subj_idx=1:size(subj_nas_session,1)
        plot(x,subj_nas_session(subj_idx,dim_idx).*1000,'o','MarkerEdgeColor',colors(subj_idx,:),'MarkerFaceColor',colors(subj_idx,:));
    end
end
for dim_idx=1:size(subj_lpa_session,2)
    x=2+(dim_idx-2)*.25;
    for subj_idx=1:size(subj_lpa_session,1)
        plot(x,subj_lpa_session(subj_idx,dim_idx).*1000,'o','MarkerEdgeColor',colors(subj_idx,:),'MarkerFaceColor',colors(subj_idx,:));
    end
end
for dim_idx=1:size(subj_rpa_session,2)
    x=3+(dim_idx-2)*.25;
    for subj_idx=1:size(subj_rpa_session,1)
        plot(x,subj_rpa_session(subj_idx,dim_idx).*1000,'o','MarkerEdgeColor',colors(subj_idx,:),'MarkerFaceColor',colors(subj_idx,:));
    end
end
xlim([0.5 3.5]);
ylabel('SD (mm)');
legend(subj_labels);
