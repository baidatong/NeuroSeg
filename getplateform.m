function plateFinal = getplateform(sigma_x,sigma_y,theta,p2)
%Get the plate weight matrix.the center of it is a plate. and the marginal
%of it is Gaussian.
h =generalGaussian(sigma_x,sigma_y,theta,p2);
% figure,mesh(h);axis off;
[x1 y1]=get_ellipse((p2-1)/2+1,(p2-1)/2+1,1.414*sigma_x,1.414*sigma_y,theta);
[X,Y]=meshgrid(1:size(h,2),1:size(h,1));
IN= inpolygon(X,Y,x1,y1);
plateform=zeros(size(h));
plateform(sub2ind(size(plateform),Y(IN),X(IN)))=1;
h2=plateform.*h;
h3=(1-plateform).*h;
plate_value=min(min(h2(IN)),max(h3(:)));
plateFinal=plate_value*plateform+h3;
sumP=sum(plateFinal(:));
if sumP~=0
    plateFinal=plateFinal/sum(plateFinal(:));
end
% figure,mesh(plateFinal);
% axis off;%axis equal;