function A=singleConv_jiasu_pre(y_cor,x_cor,I_norm,p)
%Extract the matrix of a given region. and then transform it to be a vetor.
%This is used to make convolution with the LoG matrix.
if p/2==round(p/2)
    if y_cor-p/2+1>=1 & p+y_cor+p/2<=size(I_norm,1) & x_cor-p/2+1>=1 & x_cor+p/2<=size(I_norm,2)
       aaa=I_norm(y_cor-p/2+1:y_cor+p/2,x_cor-p/2+1:x_cor+p/2);
    else
        I_norm3 = [I_norm(p:-1:1,:);I_norm;I_norm(end:-1:end-(p-1),:)];
        I_norm3 = [I_norm3(:,p:-1:1) I_norm3 I_norm3(:,end:-1:end-(p-1))];
        aaa=I_norm3(p+y_cor-p/2+1:p+y_cor+p/2,p+x_cor-p/2+1:p+x_cor+p/2);
    end
    
else %if it is on the border of the matrix
    if y_cor-floor(p/2)>=1 & y_cor+floor(p/2)<=size(I_norm,1) & x_cor-floor(p/2)>=1 & x_cor+floor(p/2)<=size(I_norm,2)
       aaa=I_norm(y_cor-floor(p/2):y_cor+floor(p/2),x_cor-floor(p/2):x_cor+floor(p/2));
    else
        I_norm3 = [I_norm(p:-1:1,:);I_norm;I_norm(end:-1:end-(p-1),:)];
        I_norm3 = [I_norm3(:,p:-1:1) I_norm3 I_norm3(:,end:-1:end-(p-1))];
        aaa=I_norm3(p+y_cor-floor(p/2):p+y_cor+floor(p/2),p+x_cor-floor(p/2):p+x_cor+floor(p/2));
    end 
end
A=aaa(:);