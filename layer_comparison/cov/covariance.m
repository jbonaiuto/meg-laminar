function covariance(fois)

spm('defaults', 'eeg')
%% arguments  
input = 'C:\empty_room\mg05125_SofieSpatialMemory_20170921_01.ds\mg05125_SofieSpatialMemory_20170921_01.meg4';
%% conversion
[dir file] = fileparts(input);
% conv_output = fullfile(dir,'SPM.dat');
% 
% S = [];
% S.dataset = input;
% S.outfile = conv_output;
% S.channels = 'MEG';
% S.timewin = [];
% S.blocksize = 3276800;
% S.checkboundary = 0;
% S.eventpadding = 0;
% S.saveorigheader = 0;
% S.conditionlabels = {'Undefined'};
% S.inputformat = [];
% S.mode = 'continuous';
% D = spm_eeg_convert(S);

D=spm_eeg_load(fullfile(dir,'SPM.mat'));
%% covariance

%Y =  squeeze(D(:,1:600,:))'; % time x chans
%co = corr(Y); % chans x chans
% co1 = co;
% figure(1)
% imagesc(co1)
nAverages = 320; 
co=zeros(size(D,1),size(D,1));
for i = 1:nAverages
    start = (i-1)*600+1;
    en =i*600;
    Y =  squeeze(D(:,start:en,:))';
    co= co+corr(Y);
end
figure()
mean_co=co./320;
new_mean_co=zeros(size(mean_co,1)+1,size(mean_co,1)+1);
new_mean_co(1:70,1:70)=mean_co(1:70,1:70);
new_mean_co(1:70,72:end)=mean_co(1:70,71:end);
new_mean_co(72:end,1:70)=mean_co(71:end,1:70);
new_mean_co(72:end,72:end)=mean_co(71:end,71:end);
new_mean_co(71,71)=1.0;
coFinal(1,:,:) = new_mean_co;
imagesc(squeeze(coFinal(1,:,:)));
set(gca,'clim',[-1 1]);
colorbar();
xlabel('Sensor');
ylabel('Sensor');
title('WB');

for j=1:size(fois,1)
    foi=fois(j,:);
    clear jobs
    matlabbatch={};
    batch_idx=1;

    % Highpass filter
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.D = {fullfile(dir,'SPM.mat')};
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.type = 'butterworth';
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.band = 'bandpass';
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.freq = foi;
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.dir = 'twopass';
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.order = 5;
    matlabbatch{batch_idx}.spm.meeg.preproc.filter.prefix = fullfile(dir,'f');
    spm_jobman('run',matlabbatch);

    D=spm_eeg_load(fullfile(dir,'fSPM.mat'));
    
%     Y =  squeeze(D(:,1:600,:))'; % time x chans
%     co = corr(Y); % chans x chans
    % co1 = co;
    % figure(1)
    % imagesc(co1)
    co=zeros(size(D,1),size(D,1));
    for i = 1:nAverages
        start = (i-1)*600+1;
        en =i*600;
        Y =  squeeze(D(:,start:en,:))';
        co= co+corr(Y);
    end
    figure()
    mean_co=co./320;
    new_mean_co=zeros(size(mean_co,1)+1,size(mean_co,1)+1);
    new_mean_co(1:70,1:70)=mean_co(1:70,1:70);
    new_mean_co(1:70,72:end)=mean_co(1:70,71:end);
    new_mean_co(72:end,1:70)=mean_co(71:end,1:70);
    new_mean_co(72:end,72:end)=mean_co(71:end,71:end);
    new_mean_co(71,71)=1.0;
    save(fullfile(dir,sprintf('cov_%d-%dHz.mat', foi(1),foi(2))),'new_mean_co');
    
    coFinal(j+1,:,:) = new_mean_co;
    imagesc(squeeze(coFinal(j+1,:,:)));
    set(gca,'clim',[-1 1]);
    colorbar();
    xlabel('Sensor');
    ylabel('Sensor');
    title(sprintf('%d-%dHz',foi(1),foi(2)));
end
    
