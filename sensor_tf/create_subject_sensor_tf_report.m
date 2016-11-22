function create_subject_sensor_tf_report(subj_info, varargin)

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
subj_dir=fullfile(params.data_dir,'analysis',subj_info.subj_id);
% Directory to put report
report_dir=fullfile(subj_dir,'sensor_tf');
if exist(report_dir,'dir')~=7
    mkdir(report_dir);
end

tpl = template('templates/','keep');
tpl = set(tpl,'file',{'subject'},{'subject_sensor_tf_template.tpl'});
tpl = set(tpl,'var',{'PAGETITLE'},{sprintf('%s: Sensor TF', subj_info.subj_id)});

file_path=fullfile('analysis',subj_info.subj_id,'time_freq_rtf_rcdots_Tafdf','t_test_positive.png');
tpl = set(tpl, 'var', {'DOTSPOSITIVETFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'time_freq_rtf_rcdots_Tafdf','t_test_negative.png');
tpl = set(tpl, 'var', {'DOTSNEGATIVETFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcdots_Tafdf','t_test_positive_broadband.png');
tpl = set(tpl, 'var', {'DOTSPOSITIVESFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcdots_Tafdf','t_test_negative_broadband.png');
tpl = set(tpl, 'var', {'DOTSNEGATIVESFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcdots_Tafdf','t_test_positive_alpha.png');
tpl = set(tpl, 'var', {'DOTSPOSITIVESFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcdots_Tafdf','t_test_negative_alpha.png');
tpl = set(tpl, 'var', {'DOTSNEGATIVESFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcdots_Tafdf','t_test_positive_beta.png');
tpl = set(tpl, 'var', {'DOTSPOSITIVESFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcdots_Tafdf','t_test_negative_beta.png');
tpl = set(tpl, 'var', {'DOTSNEGATIVESFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcdots_Tafdf','t_test_positive_gamma.png');
tpl = set(tpl, 'var', {'DOTSPOSITIVESFGAMMASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcdots_Tafdf','t_test_negative_gamma.png');
tpl = set(tpl, 'var', {'DOTSNEGATIVESFGAMMASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'time_freq_rtf_rcinstr_Tafdf','t_test_positive.png');
tpl = set(tpl, 'var', {'INSTRPOSITIVETFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'time_freq_rtf_rcinstr_Tafdf','t_test_negative.png');
tpl = set(tpl, 'var', {'INSTRNEGATIVETFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcinstr_Tafdf','t_test_positive_broadband.png');
tpl = set(tpl, 'var', {'INSTRPOSITIVESFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcinstr_Tafdf','t_test_negative_broadband.png');
tpl = set(tpl, 'var', {'INSTRNEGATIVESFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcinstr_Tafdf','t_test_positive_alpha.png');
tpl = set(tpl, 'var', {'INSTRPOSITIVESFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcinstr_Tafdf','t_test_negative_alpha.png');
tpl = set(tpl, 'var', {'INSTRNEGATIVESFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcinstr_Tafdf','t_test_positive_beta.png');
tpl = set(tpl, 'var', {'INSTRPOSITIVESFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcinstr_Tafdf','t_test_negative_beta.png');
tpl = set(tpl, 'var', {'INSTRNEGATIVESFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcinstr_Tafdf','t_test_positive_gamma.png');
tpl = set(tpl, 'var', {'INSTRPOSITIVESFGAMMASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcinstr_Tafdf','t_test_negative_gamma.png');
tpl = set(tpl, 'var', {'INSTRNEGATIVESFGAMMASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'time_freq_rtf_rcresp_Tafdf','t_test_positive.png');
tpl = set(tpl, 'var', {'RESPPOSITIVETFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'time_freq_rtf_rcresp_Tafdf','t_test_negative.png');
tpl = set(tpl, 'var', {'RESPNEGATIVETFSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcresp_Tafdf','t_test_positive_broadband.png');
tpl = set(tpl, 'var', {'RESPPOSITIVESFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcresp_Tafdf','t_test_negative_broadband.png');
tpl = set(tpl, 'var', {'RESPNEGATIVESFBROADSRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcresp_Tafdf','t_test_positive_alpha.png');
tpl = set(tpl, 'var', {'RESPPOSITIVESFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcresp_Tafdf','t_test_negative_alpha.png');
tpl = set(tpl, 'var', {'RESPNEGATIVESFALPHASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcresp_Tafdf','t_test_positive_beta.png');
tpl = set(tpl, 'var', {'RESPPOSITIVESFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcresp_Tafdf','t_test_negative_beta.png');
tpl = set(tpl, 'var', {'RESPNEGATIVESFBETASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcresp_Tafdf','t_test_positive_gamma.png');
tpl = set(tpl, 'var', {'RESPPOSITIVESFGAMMASRC'}, {[params.url_prefix file_path]});

file_path=fullfile('analysis',subj_info.subj_id,'scalp_freq_rtf_rcresp_Tafdf','t_test_negative_gamma.png');
tpl = set(tpl, 'var', {'RESPNEGATIVESFGAMMASRC'}, {[params.url_prefix file_path]});

tpl = parse(tpl,'OUT', {'subject'});

%- Display the content of the output of the parsing
out_file=fullfile(report_dir,'sensor_tf.html');
fid=fopen(out_file,'w');
%display(get(tpl,'OUT'))
fprintf(fid, strrep(strrep(get(tpl,'OUT'),'\','/'),'%','%%'));
fclose(fid);


