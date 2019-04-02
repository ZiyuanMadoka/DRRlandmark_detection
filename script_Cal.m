%% DRR pre-processing
addpath('H:/study2018/matlab_ext_libs/ReadData3D_version1k');
addpath('H:/study2018/matlab_ext_libs/MedicalImageProcessingToolbox/processing/IO');
addpath('H:/study2018/matlab_ext_libs/MedicalImageProcessingToolbox/class_image/');
addpath('H:/study2018/matlab_ext_libs/Bif/');
addpath('H:/study2018/matlab_ext_libs/encoder/encoder/');
addpath('H:/study2018/matlab_ext_libs/imshow3D/');
cd('H:/study2018/OncentraDRRs/landmark_detect/')

%% in Centos linux
addpath('/H_Disk/study2018/matlab_ext_libs/ReadData3D_version1k/');
addpath('/H_Disk/study2018/matlab_ext_libs/MedicalImageProcessingToolbox/processing/IO');
addpath('/H_Disk/study2018/matlab_ext_libs/MedicalImageProcessingToolbox/class_image/');
addpath('/H_Disk/study2018/matlab_ext_libs/Bif/');
addpath('/H_Disk/study2018/matlab_ext_libs/imshow3D/');
addpath('/H_Disk/study2018/matlab_ext_libs/encoder/encoder/');
addpath('/H_Disk/study2018/OncentraDRRs/landmark_detect/');
%% a window was already added to the cropped DRRs (T10-S1)
Initdir ='/U_Disk/Algemeen/ZWG/Pat_pydicom/auto_DRRs/';
fileset = dir(fullfile(Initdir,'*.mhd')); 
%fileset = dir('/media/ziyuan/datashare/study2017/OncentraDRRs/cropped_DRRs/pDRR_AMC*.mhd');
%ref_size = [100; 100];
dir_measure = "/U_Disk/Algemeen/ZWG/Pat_pydicom/auto_measures/measurement_";
for ii = 1:length(fileset)
    filename = fileset(ii).name;
    ntemp = split(filename,'.');
    nametemp = split(ntemp(1,1),'DRR');
    landmarks = Rib_detection(Initdir, filename,dir_measure,"N");
    write_to_file(strcat(dir_measure, nametemp(2),'.txt'),landmarks);
end

%%
function write_to_file(outputfile, LMs)
    fid = fopen(outputfile,'w');
    str_f = '%% measurements from DRR, size(mm) ,angle (degree), system coordinate system as DRR \r\n';
    formatSpec = ["Rib_width %4.2f \r\n","Length_T11L4 %4.2f \r\n", "Collimator_angle %d \r\n ","Th10_bottom_xyz %4.2f %4.2f %4.2f \r\n",...
        "L4_bottom_xyz %4.2f %4.2f %4.2f \r\n","L1_bottom_xyz %4.2f %4.2f %4.2f \r\n",...
        "Th12_right_cor_xyz %4.2f %4.2f %4.2f \r\n", "Th12_left_cor_xyz %4.2f %4.2f %4.2f \r\n",...
        "L2_right_cor_xyz %4.2f %4.2f %4.2f \r\n","L2_left_cor_xyz %4.2f %4.2f %4.2f \r\n"];
    fprintf(fid,str_f);
    %fprintf(fid,'%% measurements from DRR, scaled size 100 * 100 pixel,angle (degree), system coordinate system as DRR \r\n'); 
    %fprintf(fid,'%% the widht of the widest part of the rib \r\n');
    fprintf(fid,formatSpec(1),LMs.Rib_width);
    fprintf(fid,formatSpec(2),LMs.T11L4);
    %fprintf(fid,'%% considering bending, physical length, including Th11 and L5 \r\n');
    fprintf(fid,formatSpec(3),LMs.collimator);
    fprintf(fid,formatSpec(4),LMs.Th10_bottom_xyz);
    %fprintf(fid,'%% middle_bottem of TH10 \r\n');
    fprintf(fid,formatSpec(5),LMs.L4_bottom_xyz);
    fprintf(fid,formatSpec(6),LMs.L1_bottom_xyz);
    %fprintf(fid,'%% middle_bottem of L5 \r\n');
    fprintf(fid,formatSpec(7),LMs.Th12_right_cor_xyz);
    %fprintf(fid,'%% Right_middle of boundary of L5 \r\n');
    fprintf(fid,formatSpec(8),LMs.Th12_left_cor_xyz);
    %fprintf(fid,'%% Left_middle of boundary of Th12 \r\n');
    fprintf(fid,formatSpec(9),LMs.L2_right_cor_xyz);
    %fprintf(fid,'%% Right_middle of boundary of L5 \r\n');
    fprintf(fid,formatSpec(10),LMs.L2_left_cor_xyz);
    fclose(fid);
end
%%

figure;
subplot(3,2,1)
plot(v_Sig)
hold on
plot(rib_start:rib_last,v_Sig(rib_start:rib_last),'r')
title("1. Sum of signal along vertical direction V1")
hold off
subplot(3,2,3)
findpeaks(signal,'MinPeakWidth',margin_mm(2),'MinPeakProminence',max(signal)*0.3,'NPeaks',2,'MinPeakDistance',margin_mm(12)); 
title("2. Peak detection on half-central V1")
subplot(3,2,5)
findpeaks(-smooth(h_verte),'MinPeakDistance',I_height/12, 'NPeaks',9);
title("3. Peak detection on -H1 in vertebrae region")
subplot(3,2,[2,4,6])


%%

h = figure;axis([-10 10 -10 10]);
movie(h,frame);
set(gca,'xlim',[-1 3],'ylim',[-1 3]);
fileset2 = dir('pDRR_AMC*.mhd');


for jj = 1:length(fileset2)
    temp_fn = fileset2(jj).name;
    [img,ifo] = read_mhd(temp_fn);
    DRR = uint16(img.data);
    figure;imshow(transpose(DRR));
    frame(jj)=getframe();
    %%pause;
end

vidObj = VideoWriter('DRRs.avi');
open(vidObj);
writeVideo(vidObj,frame);
close(vidObj);
g = gabor(2,45);  %Apply the filters to the checkerboard image.
a = imgaborfilt(DRR,g);
b = fft(DRR); figure; imshow(b);
gca = figure; imshow(transpose(DRR));

%% read cut-off DRRs,
% play with different window settings
fileset = dir('pDRR_AMC*.mhd');
for ii=1 %:length(fileset)
    filename = fileset(ii).name;
    [img, ifo] = read_mhd(filename);
    DRR = transpose(uint16(img.data));
    [dim1,dim2] = size(DRR);
    T = imadjust(DRR,[0.1 0.7],[]);
    imshowpair(DRR,T,'montage');
    %imshow(transpose(img.data),[])
    %frame(ii)=getframe();
    pause;   
    %write_mhd(outputfile,img, 'NDims',2,'TransformMatrix','1 0 0 1','ElementType','uint16');
end
% try filter 
J = imnoise(DRR,'speckle');
T = imsharpen(tJ,'Amount',0.7,'Threshold',0.3);
tJ = transpose(DRR);
