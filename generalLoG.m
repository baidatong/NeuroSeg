function h =generalLoG(sigma_x,sigma_y,theta,p2)
% get a general Laplacian of Gaussian filter of size hsize with standard 
% deviation sigma_x and sigma_y. sigma_x and sigma_y are both positive.
siz = [(p2-1)/2 (p2-1)/2];
[X, Y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
a = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
b = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
c = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);

Z = exp( - (a*(X).^2 + 2*b*(X).*(Y) + c*(Y).^2)) ;
Z( Z<eps*max(Z(:))) = 0;
sumZ = sum(Z(:));
if sumZ ~= 0,
   Z  = Z/sumZ;
end
% now calculate Laplacian     
h1 = ((2*a*X+2*b*Y).^2-2*a).*Z;
h2 = ((2*b*X+2*c*Y).^2-2*c).*Z;

h = h1+h2;
h = h- sum(h(:))/prod([p2 p2]);
