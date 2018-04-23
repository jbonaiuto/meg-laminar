function plot_subject_fiducials(subj_info)

spm('defaults','eeg');

session_colors={'b','g','r','m','k','c'};
run_styles={'-','--','-.'};
legend_labels={};
figure();
hold all;

for session_num=1:length(subj_info.sessions)
    for run_num=1:subj_info.sessions(session_num)
        % Directory containing session data
        session_dir=fullfile('d:/pred_coding',subj_info.subj_id, sprintf('ses-0%d',session_num));
        % Directory containing run data
        run_code=sprintf('%s%s_JamesBonaiuto_%s_0%d', subj_info.subj_id, subj_info.birth_date, subj_info.scan_date{session_num}, run_num);
        run_dir=fullfile(session_dir, sprintf('%s.ds', run_code));
        % Directory to put results
        analysis_dir=fullfile('d:/pred_coding/derivatives/spm12', subj_info.subj_id, sprintf('ses-0%d',session_num));
        mkdir(analysis_dir);

        % File containing original MEG data
        meg_file_name=fullfile(run_dir, sprintf('%s.meg4', run_code));
        if exist(meg_file_name,'file')==2
            % File containing SPM formatted data
            spm_file_name=sprintf('spmeeg_%s.mat', run_code);

            if exist(fullfile(analysis_dir,spm_file_name),'file')~=2
                spm_jobman('initcfg');
                clear jobs
                batch_idx=1;
                matlabbatch={};

                % Convert to SPM format
                matlabbatch{batch_idx}.spm.meeg.convert.dataset = {meg_file_name};
                matlabbatch{batch_idx}.spm.meeg.convert.mode.continuous.readall = 1;
                matlabbatch{batch_idx}.spm.meeg.convert.channels{batch_idx}.all = 'all';
                matlabbatch{batch_idx}.spm.meeg.convert.outfile = fullfile(analysis_dir, spm_file_name);
                matlabbatch{batch_idx}.spm.meeg.convert.eventpadding = 0;
                matlabbatch{batch_idx}.spm.meeg.convert.blocksize = 3276800;
                matlabbatch{batch_idx}.spm.meeg.convert.checkboundary = 1;
                matlabbatch{batch_idx}.spm.meeg.convert.saveorigheader = 0;
                matlabbatch{batch_idx}.spm.meeg.convert.inputformat = 'autodetect';
                spm_jobman('run',matlabbatch);
            end

            load(fullfile(analysis_dir, spm_file_name));

            nas_x_chan=find(strcmp({D.channels.label},'HLC0011')==1);
            nas_y_chan=find(strcmp({D.channels.label},'HLC0012')==1);
            nas_z_chan=find(strcmp({D.channels.label},'HLC0013')==1);

            lpa_x_chan=find(strcmp({D.channels.label},'HLC0021')==1);
            lpa_y_chan=find(strcmp({D.channels.label},'HLC0022')==1);
            lpa_z_chan=find(strcmp({D.channels.label},'HLC0023')==1);

            rpa_x_chan=find(strcmp({D.channels.label},'HLC0031')==1);
            rpa_y_chan=find(strcmp({D.channels.label},'HLC0032')==1);
            rpa_z_chan=find(strcmp({D.channels.label},'HLC0033')==1);

            % Compute time steps
            dt=1/D.Fsample/60;
            t=[dt:dt:D.Nsamples*dt];

            % Find last event time and index
            end_evt_time=D.trials.events(end).time/60.0;
            end_evt_idx=min(find(t>=end_evt_time));

            % Get nasion, lpa, rpa location
            nas_x=D.data(nas_x_chan,:);
            nas_y=D.data(nas_y_chan,:);
            nas_z=D.data(nas_z_chan,:);

            lpa_x=D.data(lpa_x_chan,:);
            lpa_y=D.data(lpa_y_chan,:);
            lpa_z=D.data(lpa_z_chan,:);

            rpa_x=D.data(rpa_x_chan,:);
            rpa_y=D.data(rpa_y_chan,:);
            rpa_z=D.data(rpa_z_chan,:);

            % By default end index is the last event
            end_idx=end_evt_idx;

            % If there is a jump of at least 0.1 in z direction, use up until jump
            diff_z=abs(nas_z(2:end_idx)-nas_z(1:end_idx-1));
            if length(find(diff_z>0.01))
                end_idx=min(find(diff_z>0.01));
            end

            plot(nas_z(1:end_idx),sprintf('%s%s',session_colors{session_num},run_styles{run_num}));
            legend_labels{end+1}=sprintf('Session %d, Run %d', session_num, run_num);

            delete(fullfile(analysis_dir,sprintf('spmeeg_%s.mat', run_code)));
            delete(fullfile(analysis_dir,sprintf('spmeeg_%s.dat', run_code)));
        end
    end
end
legend(legend_labels);
