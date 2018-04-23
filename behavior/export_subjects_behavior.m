function export_subjects_behavior(subjects)

fid=fopen('d:/pred_coding/fig1_data.csv','w');
fprintf(fid, 'Subject,Session,Run,Trial,Congruence,Coherence,Correct,RT\n');

for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);
    for session_num=1:length(subj_info.sessions)
        for run_num=1:subj_info.sessions(session_num)
            load(fullfile('C:\pred_coding\scanning\', subj_info.subj_id, ...
                num2str(session_num), sprintf('data_%s_%d.mat', subj_info.subj_id,...
                run_num)));

            load(fullfile('C:\pred_coding\scanning\', subj_info.subj_id, ...
                num2str(session_num), sprintf('stim_%s_%d.mat', subj_info.subj_id,...
                run_num)));

            % Correct for left/right mismatch
            stim.trials(:,1)=1+2-stim.trials(:,1);
            stim.trials(:,4)=1+2-stim.trials(:,4);

            low_coherence=.5*stim.threshold;
            med_coherence=stim.threshold;
            high_coherence=1.5*stim.threshold;

            for trial_idx=1:size(stim.trials,1)
                cong='congruent';
                if stim.trials(trial_idx,3)==0
                    cong='incongruent';
                end
                coh='low';
                if stim.trials(trial_idx,2)==med_coherence
                    coh='med';
                elseif stim.trials(trial_idx,2)==high_coherence
                    coh='high';
                end
                correct=data.responses(trial_idx,1)==stim.trials(trial_idx,4);
                rt=data.responses(trial_idx,2).*1000;
                fprintf(fid, sprintf('%d,%d,%d,%d,%s,%s,%d,%.4f\n', subj_idx, session_num, run_num, trial_idx, cong, coh, correct, rt));
            end
        end
    end
end
fclose(fid);