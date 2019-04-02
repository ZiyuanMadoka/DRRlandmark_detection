function LMs = Rib_detection(Initdir, filename, dir_measure,scale)
    paren = @(x,varargin) x(varargin{:}); % indexing function-return results
    [img, ] = read_mhd(strcat(Initdir,filename));
    %DRR = uint16(img.data)';  
    DRR = img.data';
    temp_spacing = img.spacing; % pixel size in two dimensions
    margin_mm = @(length) round(length/temp_spacing(1));
    
    [img_start,I2] = precrop_filt(DRR,3); % preprocessing to remove the window and apply averaging filter 3*3 %imshow(I2); 
    [I_height, I_width] = size(I2);  % height and width of the image
    v_Sig = sum(I2(1:fix(I_height/2),:),1);  % vector_signal: sum along column, half column
    [rib_last,rib_width] = find_comp(v_Sig,1,0.25);  % select the rib region based on threshold, 0.3 between max and min
    % for UMCU_55, it was set to find_comp(v_Sig,1,0.25)
    rib_start = rib_last-rib_width+1;
    
    v_dev = abs(diff(v_Sig));
    central_start = round(I_width *0.25);
    central_end = round(I_width * 0.75);
    % minpeakdistance should be coorparated with pixel size of the image, it was 30*0.391 = 11.73 mm 
    signal =smooth(smooth(v_dev(central_start:central_end)));
    [peaks,loc]=findpeaks(signal,'MinPeakWidth',margin_mm(2),'MinPeakProminence',max(signal)*0.3,'NPeaks',2,'MinPeakDistance',margin_mm(12)); 
    
    ver_col_start =central_start+loc(1)-margin_mm(2);  % pysical size 5*0.391 = 1.955 mm
    ver_col_end = central_start+loc(2)+margin_mm(2); 
    figure;
    imshow(I2,[]); 
    title(filename,'Interpreter', 'none');
    %title('4. Landmarks on DRR')  % for SPIE paper
    hold on
    
    line([rib_start,rib_start],[1,I_height],'color','r','LineWidth',2);
    line([rib_last,rib_last],[1,I_height],'color','r','LineWidth',2);
    line([ver_col_start,ver_col_start],[1,I_height],'color','y','LineWidth',2,'LineStyle','-.');
    line([ver_col_end,ver_col_end],[1,I_height],'color','y','LineWidth',2,'LineStyle','-.');
    
    %hold off
    h_verte = sum(I2(:,ver_col_start-margin_mm(4):ver_col_end+margin_mm(4)),2); % -5 and +5  used to be 10 * 0.391 = 3.9 mm
    [vpeak1, vloc1]=findpeaks(-smooth(h_verte),'MinPeakDistance',I_height/12, 'NPeaks',9);  % 10 for not cut
    % minpeakdistance should be coorparated with pixel size of the image, it was 30*0.391 = 11.73 mm 
    vloc1 = vloc1(vloc1>10);
    if length(vloc1)<8  %was 8, if the image was cut   
        [vpeak1, vloc1]=findpeaks(diff(smooth(h_verte)),'MinPeakDistance',I_height/12, 'NPeaks',9);  % was 9, if the image was cut
        %[vpeak1, vloc1]=findpeaks(smooth(h_verte),'MinPeakDistance',round(12/temp_spacing(1)), 'NPeaks',10);  % was 9, if the image was cut
        % minpeakdistance should be coorparated with pixel size of the image, it was 30*0.391 = 11.73 mm 
    end
    vloc1 = vloc1(vloc1>10);
    vloc_sur = zeros(length(vloc1)-1,1);
    for i=1:length(vloc1)-1 
        %line([ver_col_start,ver_col_end],[vloc1(i),vloc1(i)]);
        vloc_sur(i) =  mean(sum(I2(vloc1(i):vloc1(i+1),ver_col_start-margin_mm(4):ver_col_start-margin_mm(2)),2));
    end
    length_x = min(length(vloc1),8);  % maximum 8 x points
    x1 =zeros(1,length_x);x2 = x1;vlocm = x1;
    loc_seg0 = [ver_col_start, ver_col_end];
    for i = 0:length_x-1
        xrange = ver_col_start-margin_mm(20):ver_col_end+margin_mm(20); % -25 and +25 used to be 50 * 0.391 = 19.55 mm
        minpeakdistance=margin_mm(16);
        if i==0 
            v_sig_seg = sum(I2(1:vloc1(i+1),:),1); % vector_signal: sum along column, half column
            vlocm(i+1) = (1+vloc1(i+1))/2;
        else
            if i > length(vloc1)-3
                xrange = ver_col_start-margin_mm(25):ver_col_end+margin_mm(25);
            end
            v_sig_seg = sum(I2(vloc1(i):vloc1(i+1),:),1); % vector_signal: sum along column, half column
            vlocm(i+1) = (vloc1(i)+vloc1(i+1))/2;
            minpeakdistance=max(margin_mm(16),loc_seg0(2)-loc_seg0(1)-margin_mm(10));
        end
        
        v_dev_seg = abs(diff(v_sig_seg)).* (v_sig_seg(1:length(v_sig_seg)-1)/max(v_sig_seg)).^3;
        minpeakprominence = max(smooth(v_dev_seg(ver_col_start-margin_mm(20):ver_col_end+margin_mm(20))))*0.15; %25 was 50 * 0.391 = 19.55 mm, was *0.2
       
        [peaks,loc_seg]=findpeaks(smooth(smooth(v_dev_seg(xrange))),'MinPeakWidth',5,'MinPeakProminence',minpeakprominence,'NPeaks',3,'MinPeakDistance',minpeakdistance);
        % 13 was 40 *0.391 = 15.64 m \
        idxleft = loc_seg< length(xrange)/2;   
        if(length(loc_seg(idxleft)) >= 1 && length(loc_seg(~idxleft)) >= 1)  
            [~,Ileft]=sort(peaks(idxleft),'descend');
            [~,Iright] = sort(peaks(~idxleft),'descend');
            loc_seg0 = [paren(loc_seg(idxleft),Ileft(1)),paren(loc_seg(~idxleft),Iright(1))];
            x1(i+1) = loc_seg0(1)+xrange(1);
            x2(i+1) = loc_seg0(2)+xrange(1);
            
        else
            x1(i+1) = nan;
            x2(i+1) = nan;
            loc_seg = loc_seg0;
        end
    end
    x1 = round(smooth(x1))';x2 = round(smooth(x2))';
    Fit_c3 = polyfit(vlocm(1:7),(x1(1:7)+x2(1:7))/2,1);
    plot(polyval(Fit_c3,vlocm),vlocm,'black','LineWidth',2);
    colli = atan(-Fit_c3(1))/pi*180;
    if round(colli) == 0
        collimator_angle = colli;
    else
        collimator_angle = colli/abs(colli)*floor(abs(colli));
    end
%---------------------second round for vlocs
    h_verte_2 = h_verte; %copy the size of h_verte
    if length(x1) < 8 
        for kk =(length(x1)+1):8
            x1(kk) = x1(end);
            x2(kk) = x2(end);
        end
    end
    x1_ex =[x1(1),x1];
    x2_ex = [x2(2),x2];
    for i = 1: length(vlocm)+1
        if i==1
            block = 1:vlocm(1);
        elseif i == length(vlocm)+1
            block = round(vlocm(end)+1):I_height;
        else
            block = round(vlocm(i-1)+1):vlocm(i);
        end
        h_verte_2(block) =sum(I2(block,x1_ex(i)-margin_mm(4):x2_ex(i)+margin_mm(4)),2); % -5 and +5  used to be 10 * 0.391 = 3.9 mm
    end
    [vpeak2, vloc2]=findpeaks(-smooth(h_verte_2),'MinPeakDistance',I_height/14, 'NPeaks',9);  % 10 for not cut %minpeakdistance??
    vloc2 = vloc2(vloc2>margin_mm(10));
    if i>1 && vloc2(end)>(I_height-margin_mm(6))
        Isize0 = size(I20);
        Isize = size(I2);
        vloc2(end) = vloc20(end)*Isize(1)/Isize0(1);
    end
    vloc2 = vloc2(vloc2<=(I_height-margin_mm(6))); 

    if i>1 && length(vloc2)<8
        Isize0 = size(I20);
        Isize = size(I2);
        vloc2 = vloc20*Isize(1)/Isize0(1);
    end
    %if ii == 34
    %    vloc2(4) = 172;
    %end
   
    vlocm2= vlocm;
    y1=zeros(size(x1));
    y2 = y1;
    for i=1:length(vloc2)-1
        if i==1
            vlocm2(i) = (1+vloc2(i))/2;
        else
            vlocm2(i) = (vloc2(i-1)+vloc2(i))/2;
        end
        Dy = (x2(i)-x1(i))/2*tan(Fit_c3(1)); % delta y to correct for bending
        y1(i)=vlocm2(i)+Dy;
        y2(i)=vlocm2(i)-Dy;
        line([x1(i),x2(i)],[vloc2(i)+Dy,vloc2(i)-Dy],'color','m','LineWidth',1.5);
        
        scatter(x1(i),y1(i),'MarkerFaceColor',[0.2 0.8 0.1],'MarkerEdgeColor','black'); 
        scatter(x2(i),y2(i),'MarkerFaceColor',[0.2 0.8 0.1],'MarkerEdgeColor','black');  
    % vloc_sur(i) =  mean(sum(I2(vloc2(i):vloc2(i+1),ver_col_start-margin_mm(4):ver_col_start-margin_mm(2)),2));
    end
    h_verte_sur1 =  sum(I2(1:round(vloc2(length(vloc2)-1)),ver_col_start-margin_mm(6):ver_col_start-margin_mm(2)),2);
    h_verte_sur2 =  sum(I2(1:round(vloc2(length(vloc2)-1)),ver_col_end+margin_mm(2):ver_col_end+margin_mm(4)),2);
   
    pause;
    vloc20 = vloc2;
    I20 = I2;
    close all;
    % calculate landmarks and covert back to physical size / scaled size
    if scale=="Y"
        scale_factor = ref_size./(img.size.*img.spacing); % scaled for reference
    else
        scale_factor = [1;1];
    end
    
    T11L4=sqrt(((polyval(Fit_c3,vlocm2(1))-polyval(Fit_c3,vlocm2(7)))*scale_factor(1))^2+((vlocm2(7)-vlocm2(1))*scale_factor(2))^2)*img.spacing(1);
    rib_width_fysical = sqrt((rib_width*tan(collimator_angle*pi/180)*scale_factor(2))^2+(rib_width*scale_factor(1))^2)*img.spacing(1);
    % structure to store landmarks
    LMs = struct;
    LMs.Rib_width = rib_width_fysical;
    LMs.T11L4 = T11L4;
    LMs.collimator = mod(-collimator_angle,360);
    LMs.Th10_bottom_xyz = coordi_CT(polyval(Fit_c3,vloc2(1)),vloc2(1),img,img_start);
    LMs.L4_bottom_xyz = coordi_CT(polyval(Fit_c3,vloc2(7)),vloc2(7),img,img_start);
    LMs.L1_bottom_xyz = coordi_CT(polyval(Fit_c3,vloc2(4)),vloc2(4),img,img_start);
    LMs.Th12_right_cor_xyz = coordi_CT(x1(3),y1(3),img,img_start);
    LMs.Th12_left_cor_xyz = coordi_CT(x2(3),y2(3),img,img_start);
    LMs.L2_right_cor_xyz = coordi_CT(x1(5),y1(5),img,img_start);
    LMs.L2_left_cor_xyz = coordi_CT(x2(5),y2(5),img,img_start);

end


