function max_chan=find_max_erp_channel(D, woi, erp_type, varargin)

% Parse inputs
defaults = struct('possible_channels',{});  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

num_trials=size(D,3);
good_trials=setdiff([1:num_trials],D.badtrials);
meg_ch_idx=find(strcmp(D.chantype,'MEGGRAD')==1);
if length(params.possible_channels)>0
    meg_ch_idx=[];
    for i=1:length(params.possible_channels)
        idx=find(strcmp(D.chanlabels,params.possible_channels{i})==1);
        meg_ch_idx=[meg_ch_idx idx];
    end
end

mean_D=squeeze(mean(D(:,:,good_trials),3));
ind=intersect(find(D.time>=woi(1)*0.001),find(D.time<=woi(2)*0.001));
%woi_width=woi(2)-woi(1);
%wide_ind=setdiff(intersect(find(D.time>=(woi(1)-.05)),find(D.time<=(woi(2)+.05))),ind);

if strcmp(erp_type,'P')
    %local_max_ch=[];
    %for i=1:length(meg_ch_idx)
    %    if max(mean_D(meg_ch_idx(i),ind))>max(mean_D(meg_ch_idx(i),wide_ind))
    %        local_max_ch=[local_max_ch meg_ch_idx(i)];
    %    end
    %end
    ch_max=max(mean_D(:,ind),[],2);
    %max_chan=find(ch_max==max(ch_max(intersect(meg_ch_idx,local_max_ch))));
    max_chan=find(ch_max==max(ch_max(meg_ch_idx)));
elseif strcmp(erp_type,'N')
    %local_min_ch=[];
    %for i=1:length(meg_ch_idx)
    %    if min(mean_D(meg_ch_idx(i),ind))<min(mean_D(meg_ch_idx(i),wide_ind))
    %        local_min_ch=[local_min_ch meg_ch_idx(i)];
    %    end
    %end
    ch_min=min(mean_D(:,ind),[],2);
    %max_chan=find(ch_min==min(ch_min(intersect(meg_ch_idx,local_min_ch))));
    max_chan=find(ch_min==min(ch_min(meg_ch_idx)));
end


