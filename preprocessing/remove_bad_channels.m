function channels_removed=remove_bad_channels(subj_info, session_num, data_file_name)
% function remove_bad_channels(subj_info, session_num, data_file_name)
% Remove bad channels based on scan date
% INPUT: 
%   subject info structure
%   session_num: session number
%   data_file_name: data file to remove channels from
% ---------------------------
% JJB (j.bonaiuto@ucl.ac.uk) Jul 2016
% 

spm('defaults', 'EEG');

load(data_file_name);
channels_removed={};
if datenum(subj_info.scan_date{session_num},'yyyymmdd')>=datenum('20050101','yyyymmdd')
    ch_idx=find(strcmp({D.channels.label},'MLO42')==1);
    if length(ch_idx)>0
        D.channels(ch_idx).bad=1;
        channels_removed{end+1}='MLO42';
    end
end

if datenum(subj_info.scan_date{session_num},'yyyymmdd')>=datenum('20160101','yyyymmdd')
    ch_idx=find(strcmp({D.channels.label},'MRC12')==1);
    if length(ch_idx)>0
        D.channels(ch_idx).bad=1;
        channels_removed{end+1}='MRC12';
    end
end

if datenum(subj_info.scan_date{session_num},'yyyymmdd')>=datenum('20160101','yyyymmdd')
    ch_idx=find(strcmp({D.channels.label},'MRP53')==1);
    if length(ch_idx)>0
        D.channels(ch_idx).bad=1;
        channels_removed{end+1}='MRP53';
    end
end

save(data_file_name,'D');



