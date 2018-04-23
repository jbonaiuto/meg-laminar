function compare_bias(subjects, contrasts)

out_file=fopen('C:\Users\jbonai\Dropbox\meg\pred_coding\plots\layer_comparison\bias.csv','w');
fprintf(out_file,'Subject,ROI,Contrast,FOI,Type,Z,PartialZ,Measure\n');

roi_types={'whole_brain','func','anat'};
for r_idx=1:length(roi_types)
    roi_type=roi_types{r_idx};
    disp(roi_type)
    lfn_diff_z=[];
    lfn_diff_partial_z=[];
    depth_diff_z=[];
    depth_diff_partial_z=[];

    for c_idx=1:length(contrasts)
        contrast=contrasts(c_idx);
        fid=fopen(fullfile('C:\Users\jbonai\Dropbox\meg\pred_coding\plots\layer_comparison',contrast.comparison_name,'bias',roi_type,'stats.txt'),'r');
        l_idx=1;
        tline=fgetl(fid);
        while l_idx<25
            tline=fgetl(fid);
            l_idx=l_idx+1;
        end

        s_idx=1;
        while s_idx<=length(subjects)
            tline=fgetl(fid);
            x=sscanf(tline,'LFN diff - Pial/WM diff, r=%f, p=%f');
            z=0.5*log((1+x(1))/(1-x(1)));
            lfn_diff_z(end+1,:)=[z c_idx s_idx];
            tline=fgetl(fid);
            x=sscanf(tline,'LFN diff - Pial/WM diff, partial r=%f, p=%f');
            partial_z=0.5*log((1+x(1))/(1-x(1)));
            lfn_diff_partial_z(end+1,:)=[partial_z c_idx s_idx];
            if contrast.foi(1)<=30
                band='low';
            else
                band='high';
            end            
            fprintf(out_file,'%d,%s,%s,%s,%s,%0.4f,%0.4f,LFNDiff\n', s_idx, roi_type, contrast.comparison_name, band, contrast.region, z, partial_z);
            s_idx=s_idx+1;
        end

        tline=fgetl(fid);
        tline=fgetl(fid);

        s_idx=1;
        while s_idx<=length(subjects)
            tline=fgetl(fid);
            x=sscanf(tline,'Depth diff - Pial/WM diff, r=%f, p=%f');
            z=0.5*log((1+x(1))/(1-x(1)));
            depth_diff_z(end+1,:)=[z c_idx s_idx];
            tline=fgetl(fid);
            x=sscanf(tline,'Depth diff - Pial/WM diff, partial r=%f, p=%f');
            partial_z=0.5*log((1+x(1))/(1-x(1)));
            depth_diff_partial_z(end+1,:)=[partial_z c_idx s_idx];
            if contrast.foi(1)<=30
                band='low';
            else
                band='high';
            end            
            fprintf(out_file,'%d,%s,%s,%s,%s,%0.4f,%0.4f,DepthDiff\n', s_idx, roi_type, contrast.comparison_name, band, contrast.region, z, partial_z);
            s_idx=s_idx+1;
        end

        fclose(fid);
    end
    disp('LFN diff')
    RMAOV1(lfn_diff_z,0.05)
    disp('LFN diff partial')
    RMAOV1(lfn_diff_partial_z,0.05)
    disp('Depth diff')
    RMAOV1(depth_diff_z,0.05)
    disp('Depth diff partial')
    RMAOV1(depth_diff_partial_z,0.05)
end
fclose(out_file);