function create_all_subjects_sensor_tf_report(subjects, varargin)

defaults = struct('data_dir', '/data/pred_coding', ...
    'url_prefix', 'http://fortressofjollitude.zapto.org/pred_coding/');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Add path to report generation code
addpath /home/jbonaiuto/Documents/MATLAB/m2html/
    
% Directory to put results
analysis_dir=fullfile(params.data_dir,'analysis');
% Directory to put report
report_dir=fullfile(analysis_dir,'sensor_tf');
if exist(report_dir,'dir')~=7
    mkdir(report_dir);
end

tpl = template('templates/','keep');
tpl = set(tpl,'file',{'all'},{'all_subjects_sensor_tf_template.tpl'});
tpl = set(tpl,'var',{'PAGETITLE'},{'Sensor TF'});

subj_list='<ul>';
for subj_idx=1:length(subjects)
    subj_list=sprintf('%s<li><a href="%sanalysis/%s/sensor_tf/sensor_tf.html">%s</a></li>',subj_list,params.url_prefix,subjects(subj_idx).subj_id,subjects(subj_idx).subj_id);
end
subj_list=sprintf('%s</ul>',subj_list);
tpl = set(tpl, 'var', {'SUBJLIST'}, {subj_list});

file_path=fullfile('analysis','time_freq_rtf_rcdots_Tafdf','f_test.png');
tpl = set(tpl, 'var', {'DOTSFTFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','time_freq_rtf_rcdots_Tafdf_positive','t_test_positive.png');
tpl = set(tpl, 'var', {'DOTSPOSITIVETFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','time_freq_rtf_rcdots_Tafdf_negative','t_test_negative.png');
tpl = set(tpl, 'var', {'DOTSNEGATIVETFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcdots_Tafdf','f_test_broadband.png');
tpl = set(tpl, 'var', {'DOTSFSFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcdots_Tafdf_positive','t_test_positive_broadband.png');
tpl = set(tpl, 'var', {'DOTSPOSITIVESFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcdots_Tafdf_negative','t_test_negative_broadband.png');
tpl = set(tpl, 'var', {'DOTSNEGATIVESFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcdots_Tafdf','f_test_alpha.png');
tpl = set(tpl, 'var', {'DOTSFSFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcdots_Tafdf_positive','t_test_positive_alpha.png');
tpl = set(tpl, 'var', {'DOTSPOSITIVESFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcdots_Tafdf_negative','t_test_negative_alpha.png');
tpl = set(tpl, 'var', {'DOTSNEGATIVESFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcdots_Tafdf','f_test_beta.png');
tpl = set(tpl, 'var', {'DOTSFSFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcdots_Tafdf_positive','t_test_positive_beta.png');
tpl = set(tpl, 'var', {'DOTSPOSITIVESFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcdots_Tafdf_negative','t_test_negative_beta.png');
tpl = set(tpl, 'var', {'DOTSNEGATIVESFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcdots_Tafdf','f_test_gamma.png');
tpl = set(tpl, 'var', {'DOTSFSFGAMMASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcdots_Tafdf_positive','t_test_positive_gamma.png');
tpl = set(tpl, 'var', {'DOTSPOSITIVESFGAMMASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcdots_Tafdf_negative','t_test_negative_gamma.png');
tpl = set(tpl, 'var', {'DOTSNEGATIVESFGAMMASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','time_freq_rtf_rcinstr_Tafdf','f_test.png');
tpl = set(tpl, 'var', {'INSTRFTFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','time_freq_rtf_rcinstr_Tafdf_positive','t_test_positive.png');
tpl = set(tpl, 'var', {'INSTRPOSITIVETFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','time_freq_rtf_rcinstr_Tafdf_negative','t_test_negative.png');
tpl = set(tpl, 'var', {'INSTRNEGATIVETFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcinstr_Tafdf','f_test_broadband.png');
tpl = set(tpl, 'var', {'INSTRFSFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcinstr_Tafdf_positive','t_test_positive_broadband.png');
tpl = set(tpl, 'var', {'INSTRPOSITIVESFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcinstr_Tafdf_negative','t_test_negative_broadband.png');
tpl = set(tpl, 'var', {'INSTRNEGATIVESFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcinstr_Tafdf','f_test_alpha.png');
tpl = set(tpl, 'var', {'INSTRFSFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcinstr_Tafdf_positive','t_test_positive_alpha.png');
tpl = set(tpl, 'var', {'INSTRPOSITIVESFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcinstr_Tafdf_negative','t_test_negative_alpha.png');
tpl = set(tpl, 'var', {'INSTRNEGATIVESFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcinstr_Tafdf','f_test_beta.png');
tpl = set(tpl, 'var', {'INSTRFSFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcinstr_Tafdf_positive','t_test_positive_beta.png');
tpl = set(tpl, 'var', {'INSTRPOSITIVESFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcinstr_Tafdf_negative','t_test_negative_beta.png');
tpl = set(tpl, 'var', {'INSTRNEGATIVESFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcinstr_Tafdf','f_test_gamma.png');
tpl = set(tpl, 'var', {'INSTRFSFGAMMASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcinstr_Tafdf_positive','t_test_positive_gamma.png');
tpl = set(tpl, 'var', {'INSTRPOSITIVESFGAMMASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcinstr_Tafdf_negative','t_test_negative_gamma.png');
tpl = set(tpl, 'var', {'INSTRNEGATIVESFGAMMASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','time_freq_rtf_rcresp_Tafdf','f_test.png');
tpl = set(tpl, 'var', {'RESPFTFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','time_freq_rtf_rcresp_Tafdf_positive','t_test_positive.png');
tpl = set(tpl, 'var', {'RESPPOSITIVETFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','time_freq_rtf_rcresp_Tafdf_negative','t_test_negative.png');
tpl = set(tpl, 'var', {'RESPNEGATIVETFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcresp_Tafdf','f_test_broadband.png');
tpl = set(tpl, 'var', {'RESPFSFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcresp_Tafdf_positive','t_test_positive_broadband.png');
tpl = set(tpl, 'var', {'RESPPOSITIVESFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcresp_Tafdf_negative','t_test_negative_broadband.png');
tpl = set(tpl, 'var', {'RESPNEGATIVESFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcresp_Tafdf','f_test_alpha.png');
tpl = set(tpl, 'var', {'RESPFSFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcresp_Tafdf_positive','t_test_positive_alpha.png');
tpl = set(tpl, 'var', {'RESPPOSITIVESFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcresp_Tafdf_negative','t_test_negative_alpha.png');
tpl = set(tpl, 'var', {'RESPNEGATIVESFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcresp_Tafdf','f_test_beta.png');
tpl = set(tpl, 'var', {'RESPFSFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcresp_Tafdf_positive','t_test_positive_beta.png');
tpl = set(tpl, 'var', {'RESPPOSITIVESFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcresp_Tafdf_negative','t_test_negative_beta.png');
tpl = set(tpl, 'var', {'RESPNEGATIVESFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcresp_Tafdf','f_test_gamma.png');
tpl = set(tpl, 'var', {'RESPFSFGAMMASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcresp_Tafdf_positive','t_test_positive_gamma.png');
tpl = set(tpl, 'var', {'RESPPOSITIVESFGAMMASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis','scalp_freq_rtf_rcresp_Tafdf_negative','t_test_negative_gamma.png');
tpl = set(tpl, 'var', {'RESPNEGATIVESFGAMMASRC'}, {[params.url_prefix file_path]});

tpl = parse(tpl,'OUT', {'all'});

%- Display the content of the output of the parsing
out_file=fullfile(report_dir,'sensor_tf.html');
fid=fopen(out_file,'w');
%display(get(tpl,'OUT'))
fprintf(fid, strrep(strrep(get(tpl,'OUT'),'\','/'),'%','%%'));
fclose(fid);


