function ParametersOutput=fastSegmentation(I_norm_ori,new_coordinate_x_leave,new_coordinate_y_leave,absigma_leave,MinArea,Ntheta,min_sigma)
%input : 
%         I_norm_ori              :the normalized image data
%         new_coordinate_x_leave  :the x-coordinate of the seeds
%         new_coordinate_y_leave  :the y-coordinate of the seeds
%         absigma_leave           :the sigma_x,sigma_y,and the theta of the
%                                  ellipse
%         MinArea                 :the minimum threshold of the area
%         min_sigma               :the minimum of the sigma input
%         Ntheta                  :the total portions of pi
%output: 
%       
%       ParametersOutput.xypos    :the position of the cell centroids
%       ParametersOutput.Pixels   :the position of the cells
%       ParametersOutput.xy_all   :the position of the boundaries of the 
%                                  cells
medianFilter_size=5;
try
    medianFilter = vision.MedianFilter([medianFilter_size medianFilter_size]);
    I_norm = step(medianFilter2,I_norm_ori);
catch
    I_norm = medfilt2(I_norm_ori,[medianFilter_size medianFilter_size]);% if there isn't vision package,
                                                                        % you can try medfilt2;
end
ratio_imcropSize=3;
Mask_norm=zeros(size(I_norm));
Mask_norm_temp=zeros(size(I_norm));
dataCenter=[];
edge_location=cell(0);
cell_num=0;
cell_edge=cell(0);

hwaitbar = waitbar(0,'Doing fast segmentation,please wait...');
steps=length(new_coordinate_x_leave);

for ii = 1:length(new_coordinate_x_leave)
    cell_num=ii;
    
    waitbarstep=ii;
    waitbar(waitbarstep / steps)
    
    sigma_x=absigma_leave(ii,1);
    sigma_y=absigma_leave(ii,2);
    k_temp=(absigma_leave(ii,3));
    if k_temp<0
        k_temp=0;
    end
    theta=k_temp*pi/Ntheta;
    cell_center=[new_coordinate_y_leave(ii) new_coordinate_x_leave(ii)];
    r = round(ratio_imcropSize*sqrt(2)*sigma_x); 
    imcrop=I_norm(max(1,cell_center(1)-r):min(size(I_norm,1),...
                   cell_center(1)+r),max(1,cell_center(2)-r):...
                         min(size(I_norm,2),cell_center(2)+r));
    
    cen_incrop_up=cell_center(1)-max(1,cell_center(1)-r);
    cen_incrop_left=cell_center(2)-max(1,cell_center(2)-r);
    cen_incrop_down=min(size(I_norm,1),cell_center(1)+r)-cell_center(1);
    cen_incrop_right=min(size(I_norm,2),cell_center(2)+r)-cell_center(2);
 
    [y x] = meshgrid(1:cen_incrop_up+cen_incrop_down+1,...
            1:cen_incrop_left+cen_incrop_right+1);
    y=y';
    x=x';
    
    %Caculate the weight matrix. and w_final is the final weight matrix:
    [wh2 Gx_norm Gy_norm x_norm y_norm] = imwh2(imcrop,cell_center,...
        r,cen_incrop_up,cen_incrop_left,cen_incrop_down,cen_incrop_right);
    plateFinal=getplateform(sigma_x/1.2,sigma_y/1.2,theta,2*r+1);
    
    wh1=plateFinal(r-cen_incrop_up+1:r+cen_incrop_down+1,...
        r-cen_incrop_left+1:r+cen_incrop_right+1);
    w_final=wh1.*wh2;
 
    fig_f=w_final.*imcrop;
    fig_f=(fig_f-min(fig_f(:)))/(max(fig_f(:))-min(fig_f(:)));
    
    
    %using OTSU method to segment the weighted results:
    level = graythresh(fig_f);  
    BW = im2bw(fig_f,level);
    BW2 = imfill(BW,'holes');
    imopensize=round(min_sigma);
    se = strel('disk',imopensize);
    afterOpening_patch = imopen(BW2,se);
    
%%%%%%%%%%%%%%%%%%Optional%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    se2 = strel('disk',round(imopensize/2));
    afterCloseing_patch=imclose(BW2,se2);
    afterOpening_patch = imopen(afterCloseing_patch,se);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    mask2=zeros(size(I_norm));
    mask2(new_coordinate_y_leave(ii)-cen_incrop_up:...
        new_coordinate_y_leave(ii)+cen_incrop_down,...
        new_coordinate_x_leave(ii)-cen_incrop_left:...
        new_coordinate_x_leave(ii)+cen_incrop_right)=afterOpening_patch;
    mask3=bwlabel(mask2, 4);
  
    target_num=mask3(new_coordinate_y_leave(ii),new_coordinate_x_leave(ii));
    tempcenter=[new_coordinate_y_leave(ii);new_coordinate_x_leave(ii)];
    if target_num==0
        dataCenter=[dataCenter tempcenter];
    else
        [YData XData]=find(mask3==target_num);%cleanedMask
        mask4_temp=zeros(size(mask3));
        mask4_temp(sub2ind(size(mask3),YData,XData))=1;
        YData_all{ii}=YData;
        XData_all{ii}=XData;
        XData2=[];YData2=[];
        hcon=contourc(mask4_temp,[1 1]);
        aind=hcon(2,1);
        cell_edge{cell_num}=[hcon(2,2:end);hcon(1,2:end)];
        dataCenter=[dataCenter tempcenter];
        Mask_norm=mask2 | Mask_norm;
  
    end
end
close(hwaitbar);
%%
real_cell_No=0;
STATS=cell(0);
real_cell_No_2=0;
Area=[];ConvexHull=[];Eccentricity=[];EulerNumber=[];MajorAxisLength=[];MinorAxisLength=[];ConvexArea=[];
Area_compaer_ellipse=[];
BackGround_BlurredImage_mask=1-Mask_norm;
BackGround_BlurredImage=I_norm_ori-Mask_norm.*I_norm;
ParametersOutput_old=[];
cell_intensity_imcrop_intensity=zeros(length(cell_edge),1);

hwaitbar = waitbar(0,'Determining which contour to be kept ,please wait...');
steps=length(cell_edge);

%delete some false positives:
for ii=1:length(cell_edge)
    waitbarstep=ii;
    waitbar(waitbarstep / steps);
    
    if isempty(cell_edge{ii});
        continue;
    else
      real_cell_No=real_cell_No+1;  
   
      
      Mask_prop=zeros(size(I_norm));
      Mask_prop(sub2ind(size(I_norm),YData_all{ii}',XData_all{ii}'))=1;
      
      STATS_Area= regionprops(Mask_prop,'Area');
      STATS{real_cell_No}.Area=STATS_Area.Area;
      Area=[Area;STATS{real_cell_No}.Area];

       STATS_MajorAxisLength=regionprops(Mask_prop,'MajorAxisLength');
       STATS{real_cell_No}.MajorAxisLength= STATS_MajorAxisLength.MajorAxisLength;
       MajorAxisLength=[MajorAxisLength;STATS{real_cell_No}.MajorAxisLength];
       
       STATS_MinorAxisLength=regionprops(Mask_prop,'MinorAxisLength');
       STATS{real_cell_No}.MinorAxisLength=STATS_MinorAxisLength.MinorAxisLength;
       MinorAxisLength=[MinorAxisLength;STATS{real_cell_No}.MinorAxisLength];
       
       Area_ellipse=pi*STATS_MinorAxisLength.MinorAxisLength*STATS_MajorAxisLength.MajorAxisLength/4;
       STATS{real_cell_No}.Area_ellipse=Area_ellipse;
       STATS_ConvexArea=regionprops(Mask_prop,'ConvexArea');
       STATS{real_cell_No}.ConvexArea=STATS_ConvexArea.ConvexArea;
       ConvexArea=[ConvexArea;STATS{real_cell_No}.ConvexArea];
       
       if     (STATS{real_cell_No}.Area_ellipse-STATS{real_cell_No}.Area)/STATS{real_cell_No}.Area>0.3
               continue
       elseif  STATS{real_cell_No}.Area>=MinArea
               blur_ind=sub2ind(size(I_norm),YData_all{ii}',XData_all{ii}');
               blur_sort=sort(I_norm_ori(blur_ind));
               cell_intensity=mean(blur_sort(round(length(blur_sort)/5):round(4*length(blur_sort)/5)));
              
               sigma_x=absigma_leave(ii,1);
               sigma_y=absigma_leave(ii,2);
               k_temp=(absigma_leave(ii,3));
               if k_temp<0
                  k_temp=0;
               end
               theta=k_temp*pi/Ntheta;
               r = round(ratio_imcropSize*sqrt(2)/2*sigma_x);
               BackGround_BlurredImage_temp = BackGround_BlurredImage...
                   (max(1,round(dataCenter(1,ii)-r)):min(size(I_norm,1),...
                   round(dataCenter(1,ii)+r)),max(1,...
                   round(dataCenter(2,ii)-r)):...
                   min(size(I_norm,2),round(dataCenter(2,ii)+r)));
               
               BackGround_BlurredImage_mask_temp = BackGround_BlurredImage_mask;
               imcrop_ind=find(BackGround_BlurredImage_mask_temp(max(1,round(dataCenter(1,ii)-r))...
                :min(size(I_norm,1),round(dataCenter(1,ii)+r)),max(1,round(dataCenter(2,ii)-r))...
                :min(size(I_norm,2),round(dataCenter(2,ii)+r)))==1);
               
              imcrop = BackGround_BlurredImage_temp(imcrop_ind);
               
              imcrop_sort=sort(imcrop(:));
              imcrop_intensity = mean(imcrop_sort(round(...
                   length(imcrop_sort)/4):round(3*length(imcrop_sort)/4)));
               %cell_intensity/imcrop_intensity
               cell_intensity_imcrop_intensity(ii)= cell_intensity/(imcrop_intensity+eps);

       end
    end
end
close(hwaitbar);
%%
%only save the candinate regions whose intensity larger than the background
%area
ind_in_cell=find(cell_intensity_imcrop_intensity>=1.1); 
ParametersOutput_old.xy_all=[];
for ii=1:length(ind_in_cell)
ParametersOutput_old.Pixels{1,ii}=[XData_all{ind_in_cell(ii)}';YData_all{ind_in_cell(ii)}'];
ParametersOutput_old.xy_all=[ParametersOutput_old.xy_all,[dataCenter(2,ind_in_cell(ii));dataCenter(1,ind_in_cell(ii))]];
ParametersOutput_old.xypos{1,ii}=[cell_edge{ind_in_cell(ii)}(2,:);cell_edge{ind_in_cell(ii)}(1,:)];
end
Para_ind_left1=1+zeros(length( ParametersOutput_old.xypos),1);
ParametersOutput_old1 = cell(0);
ParametersOutput_old1.xypos =  ParametersOutput_old.xypos(find(Para_ind_left1==1));
ParametersOutput_old1.Pixels = ParametersOutput_old.Pixels(find(Para_ind_left1==1));
ParametersOutput_old1.xy_all = ParametersOutput_old.xy_all(:,find(Para_ind_left1==1));
%%
%delete some cells which are overlapped
mask_all=preprocess(ParametersOutput_old1,I_norm);
ind_final=1:size(mask_all,3);
ind_final_leave=ind_final;
ind_left=[];
hwaitbar = waitbar(0,'Please wait...');

while length(ind_final)>=2
    delet_se_ind=find(mask_all(:,:,ind_final(1))==1);
    ind_com=[];
    area_com=[];
    area_com_commom=[];
    ii_ind=[];
    for ii=2:length(ind_final)            %find the overlapped area
        mask_fi_temp=mask_all(:,:,ind_final(ii));
        com_se_ind=find(mask_fi_temp(delet_se_ind)==1);
        if length(com_se_ind)>=1
            ind_com=[ind_com ind_final(ii)];
            area_com=[area_com length(find(mask_fi_temp==1))];
            area_com_commom=[area_com_commom length(com_se_ind)];
            ii_ind=[ii_ind ii];
        end
    end
    if  length(ind_com)==0 
        ind_left=[ind_left ind_final(1)];
        ind_final(1)=[];
    elseif  max(area_com_commom./length(delet_se_ind))<0.5 & max(area_com_commom./area_com) <0.5
        ind_left=[ind_left ind_final(1)];
        ind_final(1)=[]; 
    else
        if max(area_com_commom./length(delet_se_ind))>=0.5 
           ind_final(1)=[]; 
        else
            if max(area_com_commom./area_com) >=0.5                
  
              areaRatio=area_com_commom./area_com;
              area_ind_delet=find(areaRatio>=0.5);
        
              ind_final(ii_ind(area_ind_delet))=[];
              ind_left=[ind_left ind_final(1)];
              ind_final(1)=[];
            end
        end
    end
    
end
ind_left=[ind_left ind_final(1)];

ParametersOutput.Pixels = ParametersOutput_old1.Pixels(ind_left);
ParametersOutput.xy_all = ParametersOutput_old1.xy_all(:,ind_left);
ParametersOutput.xypos = ParametersOutput_old1.xypos(ind_left);
close(hwaitbar);

