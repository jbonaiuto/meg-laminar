function estimate_patch_size(subjects)

subj_dists=[];
for s_idx=1:length(subjects)
    subj_info=subjects(s_idx);
    
    D=spm_eeg_load(fullfile('D:\meg_laminar\derivatives\spm12',subj_info.subj_id,sprintf('ses-%02d', subj_info.sessions(1)),...
        'grey_coreg\EBB\p0.4\resp\f15_30',sprintf('br%s_%d.mat',subj_info.subj_id,good_sessions(1))));
    %rms=sqrt(sum(D.inv{1}.inverse.M.^2,2));
%     rms=sqrt(sum(D.inv{1}.inverse.L.^2,1));
    
    surf=gifti(fullfile('D:\meg_laminar\derivatives\freesurfer',subj_info.subj_id,...
        'ds_white.hires.deformed-ds_pial.hires.deformed.surf.gii'));

    vert  = D.inv{1}.mesh.tess_mni.vert;
    face  = D.inv{1}.mesh.tess_mni.face;
    A     = spm_mesh_distmtx(struct('vertices',vert,'faces',face),0);

    Nd=size(surf.vertices,1);

    GL    = A - spdiags(sum(A,2),0,Nd,Nd);
    GL    = GL*0.4/2;
    Qi    = speye(Nd,Nd);
    QG    = sparse(Nd,Nd);

    for i = 1:8
        QG = QG + Qi;
        Qi = Qi*GL/i;
    end

    QG    = QG.*(QG > exp(-8));
    QG    = QG*QG;
        
    % Get sample of 100 random patches
    %p_idx=randperm(D.inv{1}.inverse.Np,100);
    
    dists=[];
    for i=1:D.inv{1}.inverse.Np
        vert_idx=D.inv{1}.inverse.Ip(i);
        q=QG(:,vert_idx);
        %assert(find(q==max(q))==vert_idx);
        %rel_val=q./max(q);
        %dist=sqrt(sum((surf.vertices-repmat(surf.vertices(vert_idx,:),size(surf.vertices,1),1)).^2,2));
        val_diff=abs(.5-q./max(q));
        [xs,index]=sort(val_diff);        
        close_half_vertx=index(1:10);
        
        dists(i)=mean(sqrt(sum((surf.vertices(close_half_vertx,:)-repmat(surf.vertices(vert_idx,:),10,1)).^2,2)));
    end
           
    % For each, find closest patch
    %other_vert_idx=setdiff(D.inv{1}.inverse.Ip,vert_idx);
    %closest_verts=dsearchn(surf.vertices(other_vert_idx,:),surf.vertices(vert_idx,:));
    % Distance from each sample patch to nearest patch
    %dists=sqrt(sum(surf.vertices(other_vert_idx(closest_verts),:)-surf.vertices(vert_idx,:),2).^2);

    %mean_dist=mean(dists.*0.5);
    %subj_dists(end+1)=mean_dist;
    subj_dists(end+1)=mean(dists);
end
disp(sprintf('mean dist=%.4f, SD=%.4f', mean(subj_dists), std(subj_dists)));

% D=spm_eeg_load(fullfile('D:\pred_coding\analysis',subj_info.subj_id,1\grey_coreg\EBB\p0.4\resp\f15_30\brgb_1.mat');
% surf=gifti('D:\pred_coding\surf\gb070167-synth\surf\ds_white.hires.deformed-ds_pial.hires.deformed.surf.gii');
% 
% % Get sample of 100 random patches
% vert_idx=D.inv{1}.inverse.Ip(randperm(D.inv{1}.inverse.Np,100));
% % For each, find closest patch
% other_vert_idx=setdiff(D.inv{1}.inverse.Ip,vert_idx);
% closest_verts=dsearchn(surf.vertices(other_vert_idx,:),surf.vertices(vert_idx,:));
% % Distance from each sample patch to nearest patch
% dists=sqrt(sum(surf.vertices(other_vert_idx(closest_verts),:)-surf.vertices(vert_idx,:),2).^2);
% 
% mean_dist=mean(dists)