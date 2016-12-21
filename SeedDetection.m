function [new_coordinate_x_leave,new_coordinate_y_leave,absigma_leave]=SeedDetection(I_norm_ori,min_sigma,max_sigma,bw,Ntheta)
%input : 
%         I_norm_ori :the normalized image data
%         min_sigma  :the minimum of the sigma
%         max_sigma  :the maximum of the sigma
%         bw         :the binarization of the image data
%         Ntheta     :the total portions of pi
%output: 
%         new_coordinate_x_leave:the x-coordinate of the seeds
%         new_coordinate_y_leave:the y-coordinate of the seeds
%         absigma_leave         :the sigma_x,sigma_y,and the theta of the
%                                ellipse

medianFilter=5;
medianFilter = vision.MedianFilter([medianFilter medianFilter]);
I_norm = step(medianFilter,I_norm_ori);


norm_method=0;%choose the normalization method to use.
p2=2*max_sigma*3+1;
ii=1;
sigma_3 = [];
No_h_sum=0;%

deltasigma=1;%delta sigma 
h_sum=zeros([p2 p2 Ntheta+1]);
h_sumtemp=0;No_h=0;
alpha=0.5;
theta=-1;
I_norm_mS=adapthisteq(I_norm);
for sigma_x = min_sigma:deltasigma:max_sigma
    sigma_y = sigma_x;
    if norm_method==1
        h{ii} = (1+alpha*log(sigma_x))*(1+alpha*log(sigma_y))*generalLoG(sigma_x,sigma_y,theta,p2);
    else
        h{ii} = (sigma_x)*(sigma_y)*generalLoG(sigma_x,sigma_y,theta,p2); 
    end
   
    sigma_3 = [sigma_3;[sigma_x sigma_y theta]];
    h_sumtemp = h{ii}+ h_sumtemp;
    No_h=No_h+1;
    ii = ii+1;
end

No_h_sum = No_h_sum+1;
h_sum(:,:,No_h_sum)= h_sumtemp/No_h;

hwaitbar = waitbar(0,'Caculating LoG matrix,please wait...');
steps=Ntheta;

for  thetan=1:Ntheta
    theta=(thetan-1)/Ntheta;
    No_h=0;
    h_sumtemp=0;
     
    waitbarstep=thetan;
    waitbar(waitbarstep / steps)
    
    for sigma_x =  min_sigma+deltasigma:deltasigma:max_sigma
        for sigma_y= min_sigma:deltasigma:sigma_x-deltasigma
            
            if norm_method==1
                h{ii} = (1+alpha*log(sigma_x))*(1+alpha*log(sigma_y))*generalLoG(sigma_x,sigma_y,theta*pi,p2);
            else
                h{ii} = (sigma_x)*(sigma_y)*generalLoG(sigma_x,sigma_y,theta*pi,p2); %±ê×¼»¯
            end
            h_sumtemp=h{ii}+h_sumtemp;
            ii = ii+1; 
            No_h = No_h+1;
            sigma_3 = [sigma_3;[sigma_x sigma_y thetan-1]];
        end
    end
    No_h_sum = No_h_sum+1;
    h_sum(:,:,No_h_sum)= h_sumtemp/No_h;
end
delete(hwaitbar);
%%
%Doing convolution
hwaitbar = waitbar(0,'Doing convolution,please wait...');
steps = size(h_sum,3);

for ii=1:size(h_sum,3)  
    waitbarstep=ii;
    waitbar(waitbarstep / steps)
    MotionBlur_temp = double(imfilter(I_norm_mS,h_sum(:,:,ii),'symmetric','conv'));
    MotionBlur{ii} = (MotionBlur_temp-min(MotionBlur_temp(:)))/(max(MotionBlur_temp(:))-min(MotionBlur_temp(:)));
end
delete(hwaitbar);

MotionBlur_sum=zeros(size(I_norm));
motionblur=zeros(size(I_norm))+1;
theta_Matrix=zeros(size(I_norm_ori));
s=zeros(size(I_norm));

for ii=1:length(MotionBlur)
    MotionBlur_sum=MotionBlur_sum+MotionBlur{ii};
    s1 = double((MotionBlur{ii}<=motionblur));   %Here find the minimum value
    s2 = double((MotionBlur{ii}>motionblur)); 
    s = s+s1;
    motionblur = MotionBlur{ii}.*s1 + motionblur.*s2;  
    theta_Matrix = (zeros(size(I_norm_ori))+ii).*s1+theta_Matrix.*s2;
end

M_min = imregionalmin(MotionBlur_sum);
M_min2 = imregionalmin(motionblur);
[y_col x_col]=find(M_min==1);
[y_col2 x_col2]=find(M_min2==1);
y_col_all=[y_col;y_col2];
x_col_all=[x_col;x_col2];

se = strel('disk',round(min_sigma*1.414));
afteropening = imopen(bw,se);
afterfilling = imfill(afteropening,'holes');
hdeep=0.5;   
stdthread=1;
[mask_delete00 mask_added00]= extendDelet(MotionBlur_sum,hdeep,stdthread);
M_min_all=zeros(size(I_norm_mS));
M_min_all(sub2ind(size(I_norm_mS),y_col_all,x_col_all))=1;
seeds_candinate=(afterfilling | mask_added00).*M_min_all;%
[y_col3 x_col3]=find(seeds_candinate==1);
ll_cluster=[x_col3 y_col3];

%use the Agglomerative hierarchical cluster tree to make cluster 
dotdistance=pdist(ll_cluster);
Z=linkage(dotdistance); 
thred_distance=(max_sigma+min_sigma)/sqrt(2);  %threshold of the distance
T=cluster(Z,'cutoff',thred_distance,'criterion','distance');
new_coordinate=zeros([max(T(:)),2]);
new_coordinate_y=zeros([1,max(T(:))]);
new_coordinate_x=zeros([1,max(T(:))]);
for ii=1:max(T(:))
     l=find(T==ii);
     new_coordinate_y(ii)=round(mean(y_col3(l)));
     new_coordinate_x(ii)=round(mean(x_col3(l)));
end
for ii=-1:thetan-1
    theta_all=sigma_3(:,3);
    sigma_ind=find(ii==theta_all);
    h_sigmaFind{ii+2}=sigma_3(sigma_ind,:);
    h_Find{ii+2}=h(sigma_ind);
end

absigma=zeros(length(new_coordinate_x),3);
%%
%caculate the optimum sigma
hwaitbar = waitbar(0,'Caculating the optimum sigma,please wait...');
steps=length(new_coordinate_x);
for  ii=1:length(new_coordinate_x)
     waitbarstep=ii;
     waitbar(waitbarstep / steps)
     theta_temp=theta_Matrix(new_coordinate_y(ii),new_coordinate_x(ii));
     h_Find_temp=h_Find{theta_temp};
     h_sigmaFind_temp=h_sigmaFind{theta_temp};
     product_temp=zeros(1,length(h_Find_temp));
     XX_pre = zeros(prod([p2 p2]),length(h_Find_temp));
     h_find_pre = zeros(prod([p2 p2]),length(h_Find_temp));
     XX= singleConv_jiasu_pre(new_coordinate_y(ii),new_coordinate_x(ii),I_norm,p2); 
     for ti=1:length(h_Find_temp)
          h_find_pre(:,ti)=h_Find_temp{ti}(:);
     end
     product_temp = sum(bsxfun(@times,h_find_pre,XX)) ;
     absigma_temp_ind=find(min(product_temp(:))==product_temp);
     absigma(ii,:)=h_sigmaFind_temp(absigma_temp_ind(1),:);
end
delete(hwaitbar)
%%
%delete some false positive seeds
mIM=imfilter(I_norm,fspecial('average',p2),'replicate');
mask_delet = zeros(size(I_norm));
del_ind=zeros([length(new_coordinate_x) 1]);

hwaitbar = waitbar(0,'Detetmining which seed to be saved,please wait...');
steps=length(new_coordinate_x);

for ii=1:length(new_coordinate_x)
    waitbarstep=ii;
    waitbar(waitbarstep / steps);
    mask_delet = zeros(size(I_norm));
   % imcrop_temp_center=I_norm(new_coordinate_y(ii),new_coordinate_x(ii));
    a=absigma(ii,1);
    if a>=min_sigma+1
    b=absigma(ii,2);
    k_temp=(absigma(ii,3));
    if k_temp<0
        k_temp=0;
    end
    theta = k_temp*pi/Ntheta;
    r = round(1.5*sqrt(2)*a);
    [x1 y1] = get_ellipse(new_coordinate_x(ii),new_coordinate_y(ii),1.414*a,1.414*b,theta);
    [X,Y] = meshgrid(max(new_coordinate_x(ii)-r,1):min(new_coordinate_x(ii)+r,size(I_norm,2)),max(new_coordinate_y(ii)-r,1):min(new_coordinate_y(ii)+r,size(I_norm,1)));
    IN = inpolygon(X,Y,x1,y1);
    mask_delet(sub2ind(size(mask_delet),Y(IN),X(IN)))=1;
    dilate_se=strel('disk',round(sqrt(2)*b));
    mask_delet_dilate=imdilate(mask_delet,dilate_se);
    cen=mask_delet.*I_norm;
    back_mask=mask_delet_dilate-mask_delet;
    back=back_mask.*I_norm;
    %censort=sort(cen);
    censort_pre=sort(cen(find(mask_delet==1)));
    censort_div=censort_pre(round(length(censort_pre)/5):round(4*length(censort_pre)/5));
    backsort_pre=sort(back(find(back_mask==1)));
    backsort_div=backsort_pre(round(length(backsort_pre)/5):round(4*length(backsort_pre)/5));
    %
    del_ind(ii)=mean(censort_div(:))/mean(backsort_div(:)); 
    cell_center=[new_coordinate_y(ii) new_coordinate_x(ii)];
    imcrop=I_norm(max(1,cell_center(1)-r):min(size(I_norm,1),cell_center(1)+r),max(1,cell_center(2)-r):min(size(I_norm,2),cell_center(2)+r));
    ave_jud=mean(imcrop(:));
    std_jud=std(imcrop(:));
    a1_min=0.5;
    a1_max=3;
    a2=1;   
    center=I_norm(cell_center(1),cell_center(2));
    end
end
delete(hwaitbar)
ind_leave=find(del_ind(:)>=1.1);
new_coordinate_y_leave=new_coordinate_y(ind_leave);
new_coordinate_x_leave=new_coordinate_x(ind_leave);
absigma_leave=absigma(ind_leave,:);
end