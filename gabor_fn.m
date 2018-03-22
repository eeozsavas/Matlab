function rec=gabor_fn(I,thetas,fres)

ctr=1;
[m n]=size(I);
rec_all=zeros(m,n);

% Gabor filtering in the spatial domain
%******************************
for i=1:length(thetas)
    theta=thetas(i);
    for j=1:length(fres)
        f=fres(j);
        
        lambda=1/f;
        sigma_x=1.12*lambda;
        sigma_y = sigma_x/2;
         
        % Bounding box
        nstds = 3;
        xmax = max(abs(nstds*sigma_x*cos(theta)),abs(nstds*sigma_y*sin(theta)));
        xmax = ceil(max(1,xmax));
        ymax = max(abs(nstds*sigma_x*sin(theta)),abs(nstds*sigma_y*cos(theta)));
        ymax = ceil(max(1,ymax));
        xmin = -xmax; ymin = -ymax;
        [x,y] = meshgrid(xmin:xmax,ymin:ymax);
 
        x_theta=x*cos(theta)+y*sin(theta);
        y_theta=-x*sin(theta)+y*cos(theta);
 
        gb= exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*cos(2*pi/lambda*x_theta);
        
        res=abs(conv2(I,gb,'same'));
%         figure; imshow(res,[]);    
%         res=medfilt2(res);
        
        rec_all=rec_all+res;
        ctr=ctr+1;
        clear res gb f Ix Iy
    end
end
%******************************
rec=(rec_all)./ctr;








% % Feature selection process
% %******************************
% err_all=sum(sum(rec_all.^2));
% 
% for i=1:ctr-1
%     errors(i)=compute_errors(O(:,:,i),rec_all,err_all);
% end
% 
% [sort_errors ix]=sort(errors,'descend');
% max_error=sort_errors(1);
% rec(:,:)=O(:,:,ix(1));
% % imtool(rec,[]);


% tmpI=reshape(rec,[],1);
% ort=mean(tmpI);
% sapma=std(tmpI);
% tmpI=(tmpI-ort)/sapma;
% mi=min(tmpI);
% ma=max(tmpI);
% tmpI=(tmpI-mi)/(ma-mi);
% res=reshape(tmpI,m,n);


% OO(:,:,1)=O(:,:,ix(1));
% ind=1;
% while max_error<0.95
%     ind=ind+1;
%     rec(:,:)=rec(:,:)+O(:,:,ix(ind));
%     OO(:,:,ind)=O(:,:,ix(ind)); % Selected feature images
%     temp_err=compute_errors(rec,rec_all,err_all);
%     if temp_err>=0.95
%         break;
%     end
% end
% ctr=ind;
% 
% for i=1:ctr
%     figure; imshow(OO(:,:,i),[]);
% end