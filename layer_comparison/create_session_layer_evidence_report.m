function tpl_txt=create_session_layer_evidence_report(subj_info, session_num, zero_event, foi, comparison_name, view, varargin)

% Parse inputs
defaults = struct('data_dir', '/data/pred_coding', 'url_prefix', 'http://fortressofjollitude.zapto.org/', 'inv_type', 'EBB', 'patch_size',0.4, 'surf_dir', '');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end
if length(params.surf_dir)==0
    params.surf_dir=fullfile(params.data_dir,'surf');
end

addpath C:/users/jbonai/Documents/MATLAB/m2html/

tpl = template('templates/','keep');
tpl = set(tpl,'file',{'session'},{'session_layer_evidence_template.tpl'});
tpl = set(tpl,'var',{'TITLE'},{sprintf('Session %d', session_num)});

foi_dir=fullfile('analysis', subj_info.subj_id, num2str(session_num), 'grey_coreg', params.inv_type, ['p' num2str(params.patch_size)], zero_event, ['f' num2str(foi(1)) '_' num2str(foi(2))]);
img_path=fullfile(foi_dir,sprintf('%s.png',comparison_name));
output_file=fullfile(params.data_dir,img_path);
[pial_data,white_data,pial_white_data,pial_white_diff_data]=plot_comparison(subj_info, fullfile(params.data_dir,foi_dir), comparison_name, view, 'output_file',output_file, 'subjects_dir', params.surf_dir);
img_src=[params.url_prefix img_path];
tpl = set(tpl,'var',{'IMGSRC'},{img_src});
close all;

%- Parse template 'box' and then template 'page' and return output in 'OUT'

tpl = parse(tpl,'OUT', {'session'});

tpl_txt=get(tpl,'OUT');
