function [cond_erps trial_times session_tpl_text]=session_erp(subj_info, session_num, zero_event, woi, erp_type, varargin)

defaults = struct('data_dir', '/data/pred_coding', ...
    'url_prefix', 'http://fortressofjollitude.zapto.org/', 'time_limits',[]);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

addpath C:/users/jbonai/Documents/MATLAB/m2html/

session_tpl = template('templates/','keep');
session_tpl = set(session_tpl,'file',{'session'},{'session_erp_template.tpl'});
session_tpl = set(session_tpl,'var',{'SESSIONNUM'},{num2str(session_num)});

analysis_dir=fullfile(params.data_dir, 'analysis',subj_info.subj_id, num2str(session_num));
data_file_name=fullfile(analysis_dir, ['rc' zero_event '_Tafdf' num2str(session_num) '.mat']);

% Directory to put report
report_dir=fullfile(analysis_dir,'erp');
%delete(fullfile(report_dir,'*'));
if exist(report_dir,'dir')~=7
    mkdir(report_dir);
end

conditions={'congruent - low','congruent - med','congruent - high','incongruent - low','incongruent - med','incongruent - high'};
coherence_conditions={'low','med','high'};
congruence_conditions={'congruent','incongruent'};

cond_erps=dict();

spm('defaults', 'EEG');

session_data=spm_eeg_load(data_file_name);

possible_channels={};
if strcmp(zero_event, 'resp')
    possible_channels={'MLT13','MLT14','MLF67','MLC17'};
elseif strcmp(zero_event, 'instr')
    for i=1:length(session_data.chanlabels)
        if length(strfind(session_data.chanlabels{i},'O'))
            possible_channels{end+1}=session_data.chanlabels{i};
        end
    end
    %possible_channels={'MLO21','MLO22','MLO23','MLO24','MLO32','MLO33','MRO22','MRO32','MZO01'};
end


chan_index = find_max_erp_channel(session_data,woi,erp_type);%,'possible_channels',{possible_channels});
chan_label=session_data.chanlabels{chan_index};
session_tpl = set(session_tpl,'var',{'CHAN'},{chan_label});

% All meg channel idx
meg_ch_idx = session_data.indchantype('MEG');        
% Position of each meg channel
ch_pos=session_data.coor2D(meg_ch_idx);
% Label for each meg channel
ch_labels=session_data.chanlabels(meg_ch_idx);

vals=zeros(length(meg_ch_idx),1);
for j=1:length(meg_ch_idx)
    cur_label=ch_labels{j};
    if strcmp(cur_label,chan_label)
        vals(j)=1;
    end
end

[ZI,f]=spm_eeg_plotScalpData(vals,ch_pos,ch_labels);
%set(gca,'clim',[0 1]);
img_path=fullfile('analysis',subj_info.subj_id,num2str(session_num), 'erp',[zero_event '_channel.png']);
saveas(gcf, fullfile(params.data_dir,img_path));
close(gcf);
session_tpl = set(session_tpl,'var',{'CHANMAP'},{[params.url_prefix img_path]});

% For each coherence condition
for coher_cond_idx=1:length(coherence_conditions)

    % Get trials in this condition - congruent and incongruent
    coherence_condition=coherence_conditions{coher_cond_idx};

    for cong_cond_idx=1:length(congruence_conditions)

        congruence_condition=congruence_conditions{cong_cond_idx};
        condition=[congruence_condition ' - ' coherence_condition];
        % Get trials for this condition        
        cond_trials=session_data.indtrial(condition,'GOOD');
        disp([condition ' - ' num2str(length(cond_trials)) ' trials']);

        % If this condition has any trials
        if length(cond_trials)>1
            cond_erps(coherence_condition)=[cond_erps(coherence_condition) squeeze(session_data(chan_index,:,cond_trials))];
            cond_erps(congruence_condition)=[cond_erps(congruence_condition) squeeze(session_data(chan_index,:,cond_trials))];
            cond_erps(condition)=squeeze(session_data(chan_index,:,cond_trials));
            cond_erps('all')=[cond_erps('all') squeeze(session_data(chan_index,:,cond_trials))];
        end
    end
end   

trial_times=session_data.time([],'ms');

img_path=fullfile('analysis',subj_info.subj_id,num2str(session_num), 'erp',[zero_event '_all.png']);
plot_erp(trial_times, cond_erps, {'all'}, 'time_limits', params.time_limits, 'output_file', fullfile(params.data_dir,img_path));
session_tpl = set(session_tpl,'var',{'ALLSRC'},{[params.url_prefix img_path]});

img_path=fullfile('analysis',subj_info.subj_id,num2str(session_num), 'erp',[zero_event '_coherence.png']);
plot_erp(trial_times, cond_erps, coherence_conditions, 'time_limits', params.time_limits, 'output_file', fullfile(params.data_dir,img_path));
session_tpl = set(session_tpl,'var',{'COHERENCESRC'},{[params.url_prefix img_path]});

img_path=fullfile('analysis',subj_info.subj_id,num2str(session_num), 'erp',[zero_event '_congruence.png']);
plot_erp(trial_times, cond_erps, congruence_conditions,'time_limits', params.time_limits, 'output_file', fullfile(params.data_dir,img_path));
session_tpl = set(session_tpl,'var',{'CONGRUENCESRC'},{[params.url_prefix img_path]});

img_path=fullfile('analysis',subj_info.subj_id,num2str(session_num), 'erp',[zero_event '_congruence_coherence.png']);
plot_erp(trial_times, cond_erps, conditions, 'time_limits', params.time_limits, 'output_file', fullfile(params.data_dir,img_path));
session_tpl = set(session_tpl,'var',{'CONGRUENCECOHERENCESRC'},{[params.url_prefix img_path]});

session_tpl = parse(session_tpl,'OUT', {'session'});
session_tpl_text=strrep(strrep(get(session_tpl,'OUT'),'\','/'),'%','%%');
