function [mask_delet mask_added]= extendDelet(MotionBlur_sum,hdeep,stdthread)
%Using h-minima transform to add some regions which include some 
%candinate seeds.
mask = imextendedmin(MotionBlur_sum,hdeep);
D1 = imimposemin(MotionBlur_sum,mask);
% figure,imagesc(D1);axis equal;axis off
 inf_ind=find(D1==-inf);
 mask_inf=zeros(size(D1));
 mask_inf(inf_ind)=1;
 %figure,imagesc(mask_inf);colormap(gray);axis equal;
 CC = bwconncomp(mask_inf);
 S = regionprops(CC,'Centroid');

se_size=5;
se=strel('disk',se_size);
D1_ave=zeros([1 CC.NumObjects]);
% figure,imagesc(D1)
for ii=1:CC.NumObjects
    maskcc=zeros(size(D1));
    maskcc(CC.PixelIdxList{ii})=1;
    maskcc_dilate = imdilate(maskcc,se);
    mask_ring = maskcc_dilate-maskcc;
    mask_ring_D1=D1(find(mask_ring==1));
    mask_ring_D1(find(mask_ring_D1==-inf))=[];
    D1_ave(ii)=mean(mask_ring_D1);
%     text(S(ii).Centroid(1),S(ii).Centroid(2),num2str(ii),'color','green','FontSize',10);
end
% figure,imagesc(I_norm);colormap(gray);axis off;
% axis equal;hold on;
ll_f=find(D1~=-inf);
D1_All_ave=mean(D1(ll_f));
D1_All_std=std(D1(ll_f));
S_ind_delete=find(D1_ave(:)>=D1_All_ave-stdthread*D1_All_std);
mask_delet=ones(size(MotionBlur_sum));
for ii=1:length(S_ind_delete)
    mask_delet(CC.PixelIdxList{S_ind_delete(ii)})=0;
end
S_ind_added=find(D1_ave(:)<D1_All_ave-3/2*stdthread*D1_All_std);
mask_added=zeros(size(D1));
for ii=1:length(S_ind_added)
    mask_added(CC.PixelIdxList{S_ind_added(ii)})=1;
end

se2=strel('disk',se_size);
mask_added=imdilate(mask_added,se2);



