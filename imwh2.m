function [wh2 Gx_norm Gy_norm x_norm y_norm] = imwh2(imcrop,cell_center,r,cen_incrop_up,cen_incrop_left,cen_incrop_down,cen_incrop_right)
         %caculate the product of the gradient normal and the seed normal

         [x_norm y_norm] = ImcropdirectionNorm(cell_center,r,cen_incrop_up,cen_incrop_left,cen_incrop_down,cen_incrop_right);
         h_ground_hsize=5;
         h_ground_sigma=3;
         h_ground=fspecial('gaussian',h_ground_hsize, h_ground_sigma);
         portion2=imfilter(imcrop,h_ground);
         [Gx, Gy] = imgradientxy(portion2);
         %[Gmag, Gdir] = imgradient(Gx, Gy);
         Gm=sqrt(Gx.*Gx+Gy.*Gy);
         l=find(Gm==0);
         Gx_norm=Gx./Gm;
         Gy_norm=Gy./Gm;%quiver(x,y,Gx_norm,Gy_norm)
         Gx_norm(l)=x_norm(l);
         Gy_norm(l)=y_norm(l);
         wh2=(1+(x_norm.*Gx_norm+y_norm.*Gy_norm))/2;
end

function  [x_norm y_norm] = ImcropdirectionNorm(cell_center,r,cen_incrop_up,cen_incrop_left,cen_incrop_down,cen_incrop_right)
          %caculate the seed norm    
[y x]=meshgrid(1:cen_incrop_up+cen_incrop_down+1,1:cen_incrop_left+cen_incrop_right+1);
y=y';x=x';
x_initial=cen_incrop_left+1-x;  %(r+1)-x;
y_initial=cen_incrop_up+1-y;    %(r+1)-y;
Mag=sqrt(x_initial .* x_initial+y_initial .* y_initial);
x_norm=x_initial./Mag;
x_norm(cen_incrop_up+1,cen_incrop_left+1)=0;
y_norm=y_initial./Mag;
y_norm(cen_incrop_up+1,cen_incrop_left+1)=0;
end
