function res=FuzzyRG(I,BW)

se=[0,1,0; 1,1,1; 0,1,0];
[r c]=find(BW==1);        
mean_r=round(mean(r));
mean_c=round(mean(c));

temp=zeros(size(I));
m=0;
       
for i=mean_r-10:mean_r+10
    for j=mean_c-10:mean_c+10
%         if (I(i,j)>thrDown && I(i,j)<thrUp && BW(i,j)==1)
        if (BW(i,j)==1)
            temp(i,j)=1;
            m=m+1;
        end
    end
end

vals=I(BW==1);
mea=mean(vals);
st=std(vals);

if (m==0)
    temp(mean_r,mean_c)=1;
    mea=I(mean_r,mean_c);
    st=mea/10;
end
       
FM=exp(double((-1.*((I-mea).^2)./((st^2)*2)))); 
FM(temp==1)=1;

marker=zeros(size(I));
marker(temp==1)=1.0;

res=imreconstruct(marker,FM); 
% imtool(res,[]);
thr=0.5;
res(res<thr)=0;
res(res>=thr)=1;
res=imclose(res,se);
res=imfill(res,'holes');