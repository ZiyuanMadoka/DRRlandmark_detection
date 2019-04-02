%% DRR pre-processing
addpath('/home/ziyuan/Documents/MATLAB/ReadData3D_version1k');
addpath('/home/ziyuan/Documents/MATLAB/MedicalImageProcessingToolbox/processing/IO');
addpath('/home/ziyuan/Documents/MATLAB/MedicalImageProcessingToolbox/class_image/');
addpath('/home/ziyuan/Documents/MATLAB/Bif');
cd('/media/ziyuan/datashare/study2017/OncentraDRRs/cropped_DRRs/')

%% add a window to the cropped DRRs
fileset = dir('Cropped_AMC*.mhd');
for ii= 1:length(fileset)
    filename = fileset(ii).name;
    [img, ifo] = read_mhd(filename);
    DRR = uint16(img.data);
    [dim1,dim2] = size(DRR);
    DRRnew = uint16(zeros(620,580));
    DRRnew(40:(39+dim1),40:(39+dim2))= DRR;
    img.data = uint16(DRRnew);
    img.size = size(DRRnew);
    img.origin = img.origin(1:2);
    img.spacing = img.spacing(1:2);
    img.orientation = img.orientation(1:2,1:2);
    img.index = img.index(1:2);
    img.ndimensions=2;
    img.D = img.D(1:2,1:2);
    outputfile = ['pDRR_AMC',filename(12:end)];
    imshow(transpose(img.data))
    frame(ii)=getframe();
    pause;
    write_mhd(outputfile,img, 'NDims',2,'TransformMatrix','1 0 0 1','ElementType','uint16');
end
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

