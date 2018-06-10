function export_data(subjects)

fid=fopen('C:\meg_laminar\derivatives\spm12\sensor_tf_data_rdk.csv','w');
fprintf(fid, 'Subject,Time,Frequency,Power\n');

epoch_name='dots';
zero_evt='instr';

for subj_idx=1:length(subjects)
    results=plot_subject_tfs(subjects(subj_idx),...
        epoch_name, zero_evt, 'plot', false);
    for t_idx=1:length(results.times)
        for f_idx=1:length(results.freqs)
            fprintf(fid,'%d,%.4f,%.4f,%.4f\n', subj_idx, results.times(t_idx), results.freqs(f_idx), results.session_mean_tfs(f_idx,t_idx));
        end
    end
end
fclose(fid);

fid=fopen('C:\meg_laminar\derivatives\spm12\sensor_tf_data_instructioncue.csv','w');
fprintf(fid, 'Subject,Time,Frequency,Power\n');

epoch_name='instr';
zero_evt='instr';

for subj_idx=1:length(subjects)
    results=plot_subject_tfs(subjects(subj_idx),...
        epoch_name, zero_evt, 'plot', false);
    for t_idx=1:length(results.times)
        for f_idx=1:length(results.freqs)
            fprintf(fid,'%d,%.4f,%.4f,%.4f\n', subj_idx, results.times(t_idx), results.freqs(f_idx), results.session_mean_tfs(f_idx,t_idx));
        end
    end
end
fclose(fid);

fid=fopen('C:\meg_laminar\derivatives\spm12\sensor_tf_data_resp.csv','w');
fprintf(fid, 'Subject,Time,Frequency,Power\n');

epoch_name='resp';
zero_evt='resp';

for subj_idx=1:length(subjects)
    results=plot_subject_tfs(subjects(subj_idx),...
        epoch_name, zero_evt, 'plot', false);
    for t_idx=1:length(results.times)
        for f_idx=1:length(results.freqs)
            fprintf(fid,'%d,%.4f,%.4f,%.4f\n', subj_idx, results.times(t_idx), results.freqs(f_idx), results.session_mean_tfs(f_idx,t_idx));
        end
    end
end
fclose(fid);
