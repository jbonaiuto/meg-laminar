function run_all_session(subjects, contrasts, varargin)

defaults = struct('data_dir', 'd:/pred_coding', 'inv_type', 'EBB',...
    'patch_size',0.4, 'surf_dir', '', 'iterations',10,... %'step_size',100,...
    'thresh_percentile',80,'roi_type','mean',...
    'whole_brain', false, 'plot_ext','');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end
if length(params.surf_dir)==0
    params.surf_dir=fullfile(params.data_dir,'surf');
end
    
spm('defaults','eeg');

for i=1:length(contrasts)
    contrast=contrasts(i);
    disp(contrast.comparison_name);
    
    hemisphere=contrast.hemisphere;
    region=contrast.region;
    if params.whole_brain
        hemisphere='';
        region='';
    end
     
    thresh_type='lower';
    switch contrast.comparison_name
        case 'dots_beta_erd'
            thresh_type='upper';
        case 'dots_alpha'
            thresh_type='upper';
    end

    subj_t_vals=[];
    subj_shuffled_t_vals=[];
    
    for subj_idx=1:length(subjects)
        subj_info=subjects(subj_idx);
        disp(subj_info.subj_id);
        subj_t_vals(subj_idx).t_vals=[];
        subj_t_vals(subj_idx).num_trials=[];
        subj_shuffled_t_vals(subj_idx).t_vals=[];
        subj_shuffled_t_vals(subj_idx).num_trials=[];
        
        % Get map from white matter to pial surface
        orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
        white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');
        white_inflated_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed_inflated.surf.gii');
        orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');
        pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');
        pial_inflated_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed_inflated.surf.gii');
        pial_white_map=map_pial_to_white(white_mesh, pial_mesh, 'mapType', 'link',...
            'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh);
        white_pial_map=map_white_to_pial(white_mesh, pial_mesh, 'mapType', 'link',...
            'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh);
        pial_hemisphere_map=get_hemisphere_map(pial_mesh, orig_pial_mesh);
        white_hemisphere_map=get_hemisphere_map(white_mesh, orig_white_mesh);
        pial=gifti(pial_mesh);
        wm=gifti(white_mesh);

        nvertices=length(pial_white_map);

        pial_diff=[];
        white_diff=[];               

        foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
            num2str(subj_info.sessions(1)), 'grey_coreg', params.inv_type,...
            ['p' num2str(params.patch_size)], contrast.zero_event,...
            sprintf('f%d_%d', contrast.foi(1), contrast.foi(2)));
        lfn_filename=fullfile(foi_dir, sprintf('br%s_%d.mat',subj_info.subj_id, subj_info.sessions(1)));
        D=spm_eeg_load(lfn_filename);
        pial_lfn=sqrt(sum(D.inv{1}.inverse.L(:,nvertices+1:end).^2,1))';
        pial_lfn=spm_mesh_smooth(pial, pial_lfn, 20);
        wm_lfn=sqrt(sum(D.inv{1}.inverse.L(:,1:nvertices).^2,1))';
        wm_lfn=spm_mesh_smooth(wm, wm_lfn, 20);        
        mapped_wm_lfn=wm_lfn(pial_white_map);
        
        for session_num=1:length(subj_info.sessions)
            foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
                num2str(session_num), 'grey_coreg', params.inv_type,...
                ['p' num2str(params.patch_size)], contrast.zero_event,...
                ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))]);

            session_pial_diff=gifti(fullfile(foi_dir, sprintf('pial.%s.diff.gii', contrast.comparison_name)));
            ntrials=size(session_pial_diff.cdata,2);
            pial_diff(:,end+1:end+ntrials)=session_pial_diff.cdata(:,:);

            session_white_diff=gifti(fullfile(foi_dir, sprintf('white.%s.diff.gii', contrast.comparison_name)));
            ntrials=size(session_white_diff.cdata,2);
            white_diff(:,end+1:end+ntrials)=session_white_diff.cdata(:,:);            
        end

        trial_idx=1;
        num_trials=size(pial_diff,2);
        trial_num=10;
        while trial_num<num_trials
            for iter=1:params.iterations
                iter_trials=datasample([1:size(pial_diff,2)],trial_num,'Replace',false);
                pial_diff_iter=pial_diff(:,iter_trials);
                wm_diff_iter=white_diff(:,iter_trials);

                pial_mean=mean(pial_diff_iter,2);
                wm_mean=mean(wm_diff_iter,2);

                mapped_wm_mean=wm_mean(pial_white_map);

                prethresh_pial_mask=intersect(find(~isinf(pial_mean)), find(pial_lfn>22));
                prethresh_wm_mask=intersect(find(~isinf(wm_mean)), find(wm_lfn>19));
                prethresh_mapped_wm_mask=intersect(find(~isinf(mapped_wm_mean)), find(mapped_wm_lfn>19));

                switch hemisphere
                    case 'left'
                        prethresh_pial_mask=intersect(prethresh_pial_mask, find(pial_hemisphere_map==1));
                        prethresh_wm_mask=intersect(prethresh_wm_mask, find(white_hemisphere_map==1));
                        prethresh_mapped_wm_mask=intersect(prethresh_mapped_wm_mask, find(white_hemisphere_map(pial_white_map)==1));
                    case 'right'
                        prethresh_pial_mask=intersect(prethresh_pial_mask, find(pial_hemisphere_map==2));
                        prethresh_wm_mask=intersect(prethresh_wm_mask, find(white_hemisphere_map==2));
                        prethresh_mapped_wm_mask=intersect(prethresh_mapped_wm_mask, find(white_hemisphere_map(pial_white_map)==2));
                end

                if length(region)>0
                    switch region
                        case 'visual'
                            [lh_pial_mask,lh_wm_mask]=get_anatomical_mask('lh_visual', pial_mesh, pial_inflated_mesh,...
                                white_mesh, white_inflated_mesh, subj_info.coords('lh_visual'), 75,...
                                'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh, 'recompute', false);
                            [rh_pial_mask,rh_wm_mask]=get_anatomical_mask('rh_visual', pial_mesh, pial_inflated_mesh,...
                                white_mesh, white_inflated_mesh, subj_info.coords('rh_visual'), 75,...
                                'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh, 'recompute', false);
                            prethresh_pial_mask=intersect(prethresh_pial_mask, union(lh_pial_mask,rh_pial_mask));
                            prethresh_wm_mask=intersect(prethresh_wm_mask, union(lh_wm_mask,rh_wm_mask));
                            prethresh_mapped_wm_mask=intersect(prethresh_mapped_wm_mask, white_pial_map(prethresh_wm_mask));
                        otherwise
                            [anat_pial_mask,anat_wm_mask]=get_anatomical_mask(region, pial_mesh, pial_inflated_mesh,...
                                white_mesh, white_inflated_mesh, subj_info.coords(region), 75,...
                                'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh, 'recompute', false);
                            prethresh_pial_mask=intersect(prethresh_pial_mask, anat_pial_mask);
                            prethresh_wm_mask=intersect(prethresh_wm_mask, anat_wm_mask);
                            prethresh_mapped_wm_mask=intersect(prethresh_mapped_wm_mask, white_pial_map(anat_wm_mask));
                    end
                end

                % Find vertices greater than threshold
                if strcmp(thresh_type,'lower')
                    pial_threshold=prctile(pial_mean(prethresh_pial_mask),params.thresh_percentile);
                    wm_threshold=prctile(wm_mean(prethresh_wm_mask),params.thresh_percentile);
                    % Create pial and white masks and mapped white mask
                    pial_mask=intersect(prethresh_pial_mask, find(pial_mean>pial_threshold));
                    wm_mask=intersect(prethresh_wm_mask, find(wm_mean>wm_threshold));
                    mapped_wm_mask=intersect(prethresh_mapped_wm_mask, find(mapped_wm_mean>wm_threshold));
                % Find vertices less than threshold
                else
                    pial_threshold=prctile(pial_mean(prethresh_pial_mask),100-params.thresh_percentile);
                    wm_threshold=prctile(wm_mean(prethresh_wm_mask),100-params.thresh_percentile);
                    % Create pial and white maskss and mapped white mask
                    pial_mask=intersect(prethresh_pial_mask, find(pial_mean<pial_threshold));
                    wm_mask=intersect(prethresh_wm_mask, find(wm_mean<wm_threshold));        
                    mapped_wm_mask=intersect(prethresh_mapped_wm_mask, find(mapped_wm_mean<wm_threshold));        
                end
                mask=union(pial_mask, mapped_wm_mask);

                % Get mean pial-wm in ROI
                pial_white_diff=abs(pial_diff_iter)-abs(wm_diff_iter(pial_white_map,:));
                pial_wm_roi_diff=mean(pial_white_diff(mask,:));
                % Perform ROI t-stat
                [tstat,p]=ttest_corrected(pial_wm_roi_diff','correction',25*var(pial_wm_roi_diff));
                subj_t_vals(subj_idx).t_vals(trial_idx,iter)=tstat;                
            end
            subj_t_vals(subj_idx).num_trials(trial_idx)=trial_num;
            trial_num=trial_num*2;
            trial_idx=trial_idx+1;
        end

        for x=1:10
            trial_idx=1;
            trial_num=10;
            while trial_num<num_trials
                for iter=1:params.iterations
                    foi_dir=fullfile(params.data_dir, 'analysis', subj_info.subj_id,...
                        'grey_coreg', params.inv_type,....
                        ['p' num2str(params.patch_size)], contrast.zero_event,...
                        ['f' num2str(contrast.foi(1)) '_' num2str(contrast.foi(2))],...
                        'shuffled', num2str(iter));
                    % Get mean pial-wm in ROI
                    shuffled_pial_diff=gifti(fullfile(foi_dir,['pial.' contrast.comparison_name '.diff.gii']));
                    shuffled_white_diff=gifti(fullfile(foi_dir,['white.' contrast.comparison_name '.diff.gii']));
                    pial_diff=shuffled_pial_diff.cdata(:,:);
                    white_diff=shuffled_white_diff.cdata(:,:);

                    iter_trials=datasample([1:size(pial_diff,2)],trial_num,'Replace',false);
                    pial_diff_iter=pial_diff(:,iter_trials);
                    wm_diff_iter=white_diff(:,iter_trials);

                    pial_mean=mean(pial_diff_iter,2);
                    wm_mean=mean(wm_diff_iter,2);

                    mapped_wm_mean=wm_mean(pial_white_map);

                    prethresh_pial_mask=intersect(find(~isinf(pial_mean)), find(pial_lfn>22));
                    prethresh_wm_mask=intersect(find(~isinf(wm_mean)), find(wm_lfn>19));
                    prethresh_mapped_wm_mask=intersect(find(~isinf(mapped_wm_mean)), find(mapped_wm_lfn>19));

                    switch hemisphere
                        case 'left'
                            prethresh_pial_mask=intersect(prethresh_pial_mask, find(pial_hemisphere_map==1));
                            prethresh_wm_mask=intersect(prethresh_wm_mask, find(white_hemisphere_map==1));
                            prethresh_mapped_wm_mask=intersect(prethresh_mapped_wm_mask, find(white_hemisphere_map(pial_white_map)==1));
                        case 'right'
                            prethresh_pial_mask=intersect(prethresh_pial_mask, find(pial_hemisphere_map==2));
                            prethresh_wm_mask=intersect(prethresh_wm_mask, find(white_hemisphere_map==2));
                            prethresh_mapped_wm_mask=intersect(prethresh_mapped_wm_mask, find(white_hemisphere_map(pial_white_map)==2));
                    end

                    if length(region)>0
                        switch region
                            case 'visual'
                                [lh_pial_mask,lh_wm_mask]=get_anatomical_mask('lh_visual', pial_mesh, pial_inflated_mesh,...
                                    white_mesh, white_inflated_mesh, subj_info.coords('lh_visual'), 75,...
                                    'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh, 'recompute', false);
                                [rh_pial_mask,rh_wm_mask]=get_anatomical_mask('rh_visual', pial_mesh, pial_inflated_mesh,...
                                    white_mesh, white_inflated_mesh, subj_info.coords('rh_visual'), 75,...
                                    'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh, 'recompute', false);
                                prethresh_pial_mask=intersect(prethresh_pial_mask, union(lh_pial_mask,rh_pial_mask));
                                prethresh_wm_mask=intersect(prethresh_wm_mask, union(lh_wm_mask,rh_wm_mask));
                                prethresh_mapped_wm_mask=intersect(prethresh_mapped_wm_mask, white_pial_map(prethresh_wm_mask));
                            otherwise
                                [anat_pial_mask,anat_wm_mask]=get_anatomical_mask(region, pial_mesh, pial_inflated_mesh,...
                                    white_mesh, white_inflated_mesh, subj_info.coords(region), 75,...
                                    'origPial', orig_pial_mesh, 'origWhite', orig_white_mesh, 'recompute', false);
                                prethresh_pial_mask=intersect(prethresh_pial_mask, anat_pial_mask);
                                prethresh_wm_mask=intersect(prethresh_wm_mask, anat_wm_mask);
                                prethresh_mapped_wm_mask=intersect(prethresh_mapped_wm_mask, white_pial_map(anat_wm_mask));
                        end
                    end

                    % Find vertices greater than threshold
                    if strcmp(thresh_type,'lower')
                        pial_threshold=prctile(pial_mean(prethresh_pial_mask),params.thresh_percentile);
                        wm_threshold=prctile(wm_mean(prethresh_wm_mask),params.thresh_percentile);
                        % Create pial and white masks and mapped white mask
                        pial_mask=intersect(prethresh_pial_mask, find(pial_mean>pial_threshold));
                        wm_mask=intersect(prethresh_wm_mask, find(wm_mean>wm_threshold));
                        mapped_wm_mask=intersect(prethresh_mapped_wm_mask, find(mapped_wm_mean>wm_threshold));
                    % Find vertices less than threshold
                    else
                        pial_threshold=prctile(pial_mean(prethresh_pial_mask),100-params.thresh_percentile);
                        wm_threshold=prctile(wm_mean(prethresh_wm_mask),100-params.thresh_percentile);
                        % Create pial and white maskss and mapped white mask
                        pial_mask=intersect(prethresh_pial_mask, find(pial_mean<pial_threshold));
                        wm_mask=intersect(prethresh_wm_mask, find(wm_mean<wm_threshold));        
                        mapped_wm_mask=intersect(prethresh_mapped_wm_mask, find(mapped_wm_mean<wm_threshold));        
                    end
                    mask=union(pial_mask, mapped_wm_mask);

                    % Get mean pial-wm in ROI
                    pial_white_diff=abs(pial_diff_iter)-abs(wm_diff_iter(pial_white_map,:));
                    pial_wm_roi_diff=mean(pial_white_diff(mask,:));
                    [tstat,p]=ttest_corrected(pial_wm_roi_diff','correction',25*var(pial_wm_roi_diff));
                    subj_shuffled_t_vals(subj_idx).t_vals(x,trial_idx,iter)=tstat;                
                end
                subj_shuffled_t_vals(subj_idx).num_trials(trial_idx)=trial_num;
                trial_num=trial_num*2;
                trial_idx=trial_idx+1;
            end
        end
        t_vals=squeeze(mean(subj_shuffled_t_vals(subj_idx).t_vals,3));
        if min(t_vals(:,end))>0
            min_idx=find(t_vals(:,end)==min(t_vals(:,end)));
            subj_shuffled_t_vals(subj_idx).t_vals=squeeze(subj_shuffled_t_vals(subj_idx).t_vals(min_idx,:,:));
        else
            if max(t_vals(:,end))>0
                pos_vals=t_vals(find(t_vals(:,end)>0),end);
                min_idx=find(t_vals(:,end)==min(pos_vals));
                subj_shuffled_t_vals(subj_idx).t_vals=squeeze(subj_shuffled_t_vals(subj_idx).t_vals(min_idx,:,:));
            else
                max_idx=find(t_vals(:,end)==max(t_vals(:,end)));
                subj_shuffled_t_vals(subj_idx).t_vals=squeeze(subj_shuffled_t_vals(subj_idx).t_vals(max_idx,:,:));
            end
        end
    end

    plot_dir=fullfile('C:\Users\jbonai\Dropbox\meg\pred_coding\plots\layer_comparison',...
            contrast.comparison_name, 'noise');
    mkdir(plot_dir);
    fig=figure();
    hold all;
    max_num_trials=0;
    alpha=1.0-(0.05/2);
    for subj_idx=1:length(subjects)
        t_vals=subj_t_vals(subj_idx).t_vals;
        num_trials=subj_t_vals(subj_idx).num_trials;
        max_num_trials=max([max_num_trials, max(num_trials)]);
        mean_tvals=mean(t_vals,2);
        std_tvals=std(t_vals,[],2)./sqrt(params.iterations);
        errorbar(num_trials,mean_tvals,std_tvals,'--','Color','k');                
        
        t_vals=subj_shuffled_t_vals(subj_idx).t_vals;
        num_trials=subj_shuffled_t_vals(subj_idx).num_trials;
        max_num_trials=max([max_num_trials, max(num_trials)]);
        mean_tvals=mean(t_vals,2);
        std_tvals=std(t_vals,[],2)./sqrt(params.iterations);
        errorbar(num_trials,mean_tvals,std_tvals,'--','Color',[0.5 0.5 0.5]);                
    end
    t_thresh=tinv(alpha, [2:max_num_trials+num_trials(1)]-1);        
    plot([2:max_num_trials+num_trials(1)],t_thresh,'k--');
    plot([2:max_num_trials+num_trials(1)],-t_thresh,'k--');    
    xlabel('Num trials');
    ylabel('Pial-White t');
    xlim([0 max_num_trials+num_trials(1)]);
    title(contrast.comparison_name);
    saveas(fig, fullfile(plot_dir, sprintf('%s_%s.fig', contrast.comparison_name, params.plot_ext)), 'fig');
    saveas(fig, fullfile(plot_dir, sprintf('%s_%s.png', contrast.comparison_name, params.plot_ext)), 'png');
    saveas(fig, fullfile(plot_dir, sprintf('%s_%s.eps', contrast.comparison_name, params.plot_ext)), 'eps');
end
