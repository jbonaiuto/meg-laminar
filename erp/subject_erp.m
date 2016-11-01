function [cond_erps trial_times]=subject_erp(subj_info, label, zero_event, woi, erp_type, varargin)

defaults = struct('data_dir', '/data/pred_coding', ...
    'url_prefix', 'http://fortressofjollitude.zapto.org/', 'time_limits',[]);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

addpath C:/users/jbonai/Documents/MATLAB/m2html/

subject_tpl = template('templates/','keep');
subject_tpl = set(subject_tpl,'file',{'subject'},{'subject_erp_template.tpl'});
subject_tpl = set(subject_tpl,'var',{'SUBJECT'},{subj_info.subj_id});

subject_tpl = set(subject_tpl,'block','subject','session','sessions');
subject_tpl = set(subject_tpl, 'var', {'PAGETITLE','ZEROEVENT','WOI'}, ...
    {['ERP: Subject ' subj_info.subj_id ': ' label], label, ...
    [erp_type num2str(woi(1)) '-' num2str(woi(2)) 's']});

% Directory to put report
report_dir=fullfile(params.data_dir, 'analysis',subj_info.subj_id,'erp');
%delete(fullfile(report_dir,'*'));
if exist(report_dir,'dir')~=7
    mkdir(report_dir);
end

conditions={'congruent - low','congruent - med','congruent - high',...
    'incongruent - low','incongruent - med','incongruent - high'};
coherence_conditions={'low','med','high'};
congruence_conditions={'congruent','incongruent'};

cond_erps=dict();
for session_num=1:length(subj_info.sessions)
    [session_erps trial_times session_text]=session_erp(subj_info, ...
        session_num, label, zero_event, woi, erp_type, 'data_dir', params.data_dir,...
        'url_prefix', params.url_prefix, 'time_limits', params.time_limits);
    keys=session_erps.keys;
    for i=1:length(keys)
        cond_erps(keys{i})=[cond_erps(keys{i}) mean(session_erps(keys{i}),2)];
    end
    subject_tpl = set(subject_tpl, 'var', {'SESSIONOUT'}, {session_text});
    subject_tpl = parse(subject_tpl,'sessions','session',1);
end

img_path=fullfile('analysis',subj_info.subj_id,'erp',[label '_all.png']);
plot_erp(trial_times, cond_erps, {'all'}, 'time_limits', params.time_limits, 'output_file', fullfile(params.data_dir,img_path));
subject_tpl = set(subject_tpl,'var',{'ALLSRC'},{[params.url_prefix img_path]});

img_path=fullfile('analysis',subj_info.subj_id,'erp',[label '_coherence.png']);
plot_erp(trial_times, cond_erps, coherence_conditions, 'time_limits', params.time_limits, 'output_file', fullfile(params.data_dir,img_path));
subject_tpl = set(subject_tpl,'var',{'COHERENCESRC'},{[params.url_prefix img_path]});

img_path=fullfile('analysis',subj_info.subj_id,'erp',[label '_congruence.png']);
plot_erp(trial_times, cond_erps, congruence_conditions, 'time_limits', params.time_limits, 'output_file', fullfile(params.data_dir,img_path));
subject_tpl = set(subject_tpl,'var',{'CONGRUENCESRC'},{[params.url_prefix img_path]});

img_path=fullfile('analysis',subj_info.subj_id,'erp',[label '_congruence_coherence.png']);
plot_erp(trial_times, cond_erps, conditions, 'time_limits', params.time_limits, 'output_file', fullfile(params.data_dir,img_path));
subject_tpl = set(subject_tpl,'var',{'CONGRUENCECOHERENCESRC'},{[params.url_prefix img_path]});

subject_tpl = parse(subject_tpl,'OUT', {'subject'});

%- Display the content of the output of the parsing
out_file=fullfile(report_dir,[label '_erp.html']);
fid=fopen(out_file,'w');
%display(get(tpl,'OUT'))
fprintf(fid, strrep(strrep(get(subject_tpl,'OUT'),'\','/'),'%','%%'));
fclose(fid);

