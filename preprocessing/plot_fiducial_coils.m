function plot_fiducial_coils(data_file_name, varargin)
% function plot_fiducial_coils(data_file_name, varargin)
% Plot fiducial coil locations over a run
% INPUT: 
%   data_file_name: name of the file containing the MEG data
%   output_file: name of file to save plot to (optional)
%   output_format: format to save plot image to (optional, default=png)
% ---------------------------
% JJB (j.bonaiuto@ucl.ac.uk) Jul 2016
% 

% Parse inputs
defaults = struct('output_file','','output_format','png','lim_last_event',true, 'lim_jump',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

load(data_file_name);

fiducials=[];

fiducials(1).name='nasion';
fiducials(1).x_chan=find(strcmp({D.channels.label},'HLC0011')==1);
fiducials(1).y_chan=find(strcmp({D.channels.label},'HLC0012')==1);
fiducials(1).z_chan=find(strcmp({D.channels.label},'HLC0013')==1);

fiducials(2).name='lpa';
fiducials(2).x_chan=find(strcmp({D.channels.label},'HLC0021')==1);
fiducials(2).y_chan=find(strcmp({D.channels.label},'HLC0022')==1);
fiducials(2).z_chan=find(strcmp({D.channels.label},'HLC0023')==1);

fiducials(3).name='rpa';
fiducials(3).x_chan=find(strcmp({D.channels.label},'HLC0031')==1);
fiducials(3).y_chan=find(strcmp({D.channels.label},'HLC0032')==1);
fiducials(3).z_chan=find(strcmp({D.channels.label},'HLC0033')==1);

% Time step
dt=1/D.Fsample/60;
end_time=D.Nsamples/D.Fsample/60.0;
t=[dt:dt:end_time];
end_idx=size(D.data,2);
end_evt_time=D.trials.events(end).time/60.0;

diffz=D.data(fiducials(1).z_chan,2:end_idx)-D.data(fiducials(1).z_chan,1:end_idx-1);
if length(find(diffz>0.01))
    jump_idx=min(find(diffz>0.01));
    jump_time=t(jump_idx);
else
    jump_idx=end_idx;
    jump_time=end_evt_time;
end

if params.lim_last_event
    % Time of end event
    end_time=end_evt_time;
    end_idx=min(find(t>=end_evt_time));
    t=t(1:end_idx);
end

if params.lim_jump
    end_time=jump_time;
    end_idx=jump_idx;
    t=t(1:end_idx);
end

fig=figure();
for idx=1:length(fiducials)
    subplot(3,1,idx);
    hold all;

    mean_x=mean(D.data(fiducials(idx).x_chan,1:100:end_idx));
    zero_mean_x=D.data(fiducials(idx).x_chan,1:end_idx)-mean_x;
    plot(t,zero_mean_x,'r');
    
    mean_y=mean(D.data(fiducials(idx).y_chan,1:100:end_idx));
    zero_mean_y=D.data(fiducials(idx).y_chan,1:end_idx)-mean_y;
    plot(t,zero_mean_y,'g');
    
    mean_z=mean(D.data(fiducials(idx).z_chan,1:100:end_idx));
    zero_mean_z=D.data(fiducials(idx).z_chan,1:end_idx)-mean_z;
    plot(t,zero_mean_z,'b');
    
    xlim([dt end_time]);
    min_y=min([min(D.data(fiducials(idx).x_chan,1:10:end_idx)-mean_x) min(D.data(fiducials(idx).y_chan,1:10:end_idx)-mean_y) min(D.data(fiducials(idx).z_chan,1:10:end_idx)-mean_z)]);
    max_y=max([max(D.data(fiducials(idx).x_chan,1:10:end_idx)-mean_x) max(D.data(fiducials(idx).y_chan,1:10:end_idx)-mean_y) max(D.data(fiducials(idx).z_chan,1:10:end_idx)-mean_z)]);
    plot([end_evt_time end_evt_time],[min_y max_y],'k--');
    plot([jump_time jump_time],[min_y max_y],'m--');
    
    hold off;
    legend({'x','y','z','end evt','jump'});
    xlabel('Time (min)');
    ylabel([fiducials(idx).name ' (m)']);
end

if length(params.output_file)>0
    saveas(fig, params.output_file, params.output_format);
end
