function extract_fiducials(subjects)

spm('defaults','eeg');

for subj_idx=1:length(subjects)
    subj_info=subjects(subj_idx);

    for session_num=1:length(subj_info.sessions)
        for run_num=1:subj_info.sessions(session_num)
            % Directory containing session data
            session_dir=fullfile('c:/meg_laminar', subj_info.subj_id, sprintf('ses-%02d',session_num),'meg');
            % Directory containing run data
            run_code=sprintf('%s%s_JamesBonaiuto_%s_0%d', subj_info.subj_id, subj_info.scan_date{session_num}, run_num);
            run_dir=fullfile(session_dir, sprintf('%s.ds', run_code));
            % Directory to put results
            analysis_dir=fullfile('d:/meg_laminar/derivatives/spm12', subj_info.subj_id, sprintf('ses-%02d',session_num));
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

                % Find second to last event time and index
                end_evt_time=D.trials.events(end-1).time/60.0;
                end_evt_idx=min(find(t>=end_evt_time));

                % Get nasion, lpa, rpa location
                mvmt.nas_x=D.data(nas_x_chan,1:end_evt_idx);
                mvmt.nas_y=D.data(nas_y_chan,1:end_evt_idx);
                mvmt.nas_z=D.data(nas_z_chan,1:end_evt_idx);

                mvmt.lpa_x=D.data(lpa_x_chan,1:end_evt_idx);
                mvmt.lpa_y=D.data(lpa_y_chan,1:end_evt_idx);
                mvmt.lpa_z=D.data(lpa_z_chan,1:end_evt_idx);

                mvmt.rpa_x=D.data(rpa_x_chan,1:end_evt_idx);
                mvmt.rpa_y=D.data(rpa_y_chan,1:end_evt_idx);
                mvmt.rpa_z=D.data(rpa_z_chan,1:end_evt_idx);

                save(fullfile(analysis_dir, sprintf('mvmt_%d.mat',run_num)), 'mvmt');

                delete(fullfile(analysis_dir,sprintf('spmeeg_%s.mat', run_code)));
                delete(fullfile(analysis_dir,sprintf('spmeeg_%s.dat', run_code)));
            end
        end
    end
end
