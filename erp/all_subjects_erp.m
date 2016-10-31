function all_subjects_erp(subjects, zero_event, woi, erp_type, varargin)

defaults = struct('data_dir', '/data/pred_coding', ...
    'url_prefix', 'http://fortressofjollitude.zapto.org/', 'time_limits',[]);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

addpath C:/users/jbonai/Documents/MATLAB/m2html/

all_tpl = template('templates/','keep');
all_tpl = set(all_tpl,'file',{'all'},{'all_subjects_erp_template.tpl'});

all_tpl = set(all_tpl,'block','all','subject','subjects');
all_tpl = set(all_tpl, 'var', {'PAGETITLE','ZEROEVENT','WOI'}, {['ERP: All subjects: ' zero_event], zero_event, [erp_type num2str(woi(1)) '-' num2str(woi(2)) 'ms']});

% Directory to put report
report_dir=fullfile(params.data_dir,'analysis','erp');
%delete(fullfile(report_dir,'*'));
if exist(report_dir,'dir')~=7
    mkdir(report_dir);
end

conditions={'congruent - low','congruent - med','congruent - high','incongruent - low','incongruent - med','incongruent - high'};
coherence_conditions={'low','med','high'};
congruence_conditions={'congruent','incongruent'};

cond_erps=dict();

for subject_idx=1:length(subjects)
    [subject_erps trial_times]=subject_erp(subjects(subject_idx), zero_event, woi, erp_type, 'data_dir', params.data_dir, 'url_prefix', params.url_prefix, 'time_limits', params.time_limits);
    keys=subject_erps.keys;
    for i=1:length(keys)
        cond_erps(keys{i})=[cond_erps(keys{i}) mean(subject_erps(keys{i}),2)];
    end
    subject_text=['<h2><a href="' params.url_prefix 'analysis/' subjects(subject_idx).subj_id '/erp/' zero_event '_erp.html">Subject ' subjects(subject_idx).subj_id '</a></h2><table><tr><td><strong>All</strong></td><td><strong>Congruence</strong></td><td><strong>Coherence</strong></td><td><strong>Congruence x Coherence</strong></td></tr><tr><td><img src="' params.url_prefix 'analysis/' subjects(subject_idx).subj_id '/erp/' zero_event '_all.png" width="100%"/></td><td><img src="' params.url_prefix 'analysis/' subjects(subject_idx).subj_id '/erp/' zero_event '_congruence.png" width="100%"/></td><td><img src="' params.url_prefix 'analysis/' subjects(subject_idx).subj_id '/erp/' zero_event '_coherence.png" width="100%"/></td><td><img src="' params.url_prefix 'analysis/' subjects(subject_idx).subj_id '/erp/' zero_event '_congruence_coherence.png" width="100%"/></td></tr></table>'];
    
    all_tpl = set(all_tpl, 'var', {'SUBJECTOUT'}, {subject_text});
    all_tpl = parse(all_tpl,'subjects','subject',1);
end

img_path=fullfile('analysis','erp',[zero_event '_all.png']);
plot_erp(trial_times, cond_erps, {'all'}, 'time_limits', params.time_limits, 'output_file', fullfile(params.data_dir,img_path));
all_tpl = set(all_tpl,'var',{'ALLSRC'},{[params.url_prefix img_path]});

img_path=fullfile('analysis','erp',[zero_event '_coherence.png']);
plot_erp(trial_times, cond_erps, coherence_conditions, 'time_limits', params.time_limits, 'output_file', fullfile(params.data_dir,img_path));
all_tpl = set(all_tpl,'var',{'COHERENCESRC'},{[params.url_prefix img_path]});

img_path=fullfile('analysis','erp',[zero_event '_congruence.png']);
plot_erp(trial_times, cond_erps, congruence_conditions, 'time_limits', params.time_limits, 'output_file', fullfile(params.data_dir,img_path));
all_tpl = set(all_tpl,'var',{'CONGRUENCESRC'},{[params.url_prefix img_path]});

img_path=fullfile('analysis','erp',[zero_event '_congruence_coherence.png']);
plot_erp(trial_times, cond_erps, conditions, 'time_limits', params.time_limits, 'output_file', fullfile(params.data_dir,img_path));
all_tpl = set(all_tpl,'var',{'CONGRUENCECOHERENCESRC'},{[params.url_prefix img_path]});


all_tpl = parse(all_tpl,'OUT', {'all'});

%- Display the content of the output of the parsing
out_file=fullfile(report_dir,[zero_event '_erp.html']);
fid=fopen(out_file,'w');
%display(get(tpl,'OUT'))
fprintf(fid, strrep(strrep(get(all_tpl,'OUT'),'\','/'),'%','%%'));
fclose(fid);
