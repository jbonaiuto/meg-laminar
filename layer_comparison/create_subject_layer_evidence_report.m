function tpl_txt=create_subject_layer_evidence_report(subj_info, zero_event, foi, comparison_name, view, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'url_prefix', 'http://fortressofjollitude.zapto.org/pred_coding/', 'inv_type', 'EBB', 'patch_size',0.4, 'surf_dir', '');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end
if length(params.surf_dir)==0
    params.surf_dir=fullfile(params.data_dir,'surf');
end

addpath /home/jbonaiuto/Documents/MATLAB/m2html/

tpl = template('templates/','keep');
tpl = set(tpl,'file',{'subject'},{'subject_layer_evidence_template.tpl'});
tpl = set(tpl,'block','subject','session','sessions');
tpl = set(tpl,'var',{'TITLE'},{sprintf('Subject %s, %s', subj_info.subj_id, comparison_name)});

for i=1:length(subj_info.sessions)
    session_txt=create_session_layer_evidence_report(subj_info, i, zero_event, foi, comparison_name, view, 'data_dir', params.data_dir, 'url_prefix', params.url_prefix, 'inv_type', params.inv_type, 'patch_size', params.patch_size, 'surf_dir', params.surf_dir);
    tpl = set(tpl, 'var', {'SESSIONOUT'}, {session_txt});
    tpl = parse(tpl,'sessions','session',1);

end

tpl = parse(tpl,'OUT', {'subject'});
%- Display the content of the output of the parsing
report_dir=fullfile(params.data_dir,'analysis',subj_info.subj_id,'grey_coreg', params.inv_type, ['p' num2str(params.patch_size)], zero_event, ['f' num2str(foi(1)) '_' num2str(foi(2))]);
mkdir(report_dir);
out_file=fullfile(report_dir,sprintf('layer_evidence_%s.html',comparison_name));
fid=fopen(out_file,'w');
%display(get(tpl,'OUT'))
tpl_txt=strrep(strrep(get(tpl,'OUT'),'\','/'),'%','%%');
fprintf(fid, tpl_txt);
fclose(fid);

