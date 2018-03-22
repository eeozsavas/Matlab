% spinal canal segmentasyonu icin kullanilan fonksiyon

function ind = findFirstClosedSp(sliceNumber)

Z1=zeros(512,512,sliceNumber);
inVolume=zeros(512,512,sliceNumber);

ftitle='HU.mat';
load(ftitle);

ftitle='results.mat';
load(ftitle);

for i=sliceNumber:-1:1
    XX(:,:)=inVolume(:,:,i);
    Z(:,:)=Z1(:,:,i);
    [s1 s2]=size(Z);
    Sp=zeros(s1,s2);
    
    %patient'ý içeren pencere
    %***********************************************************
    [r c]=find(Z~=0);
    minc=min(c);
    maxc=max(c);
    minr=min(r);
    maxr=max(r);
    clear r c
    centerc=round((minc+maxc)/2);
    centerr=round((minr+maxr)/2);
    %***********************************************************
    
    Y=XX;
    Y(Z==0)=-10000; %sadece patient kýsýmlarýný alýyor

    [s1 s2]=size(Y);        
    Bone=zeros(s1,s2);
    Bone(Z==4)=1;
    clear XX Z

    %patient kýsmýnýn ortasýndan geçen çizgiye deðen kemikler vertebra
    %olarak iþaretleniyor
    %***********************************************************
    [L num]=bwlabel(Bone);
    Mask=zeros(s1,s2);
    Mask(centerr-50:maxr-10,centerc-10:centerc+10)=1;
    V=zeros(s1,s2);

    for ii=1:num
        A=zeros(s1,s2);
        A(L==ii)=1;
        if (length(nonzeros(A))>50)
            temp=A.*Mask;
            if (length(nonzeros(temp))>0)
                V(L==ii)=1;
            end
        end
        clear A temp
    end

    E=zeros(s1,s2);
    [px,py] = gradient(V);
    mag=(px.^2+py.^2).^0.5;
    for ix=1:s1
        for jy=1:s2
            if (mag(ix,jy)~=0)
                E(ix,jy)=1;
            end
        end
    end
    V(E==1)=1;
    V(:,1:centerc-50)=0;
    V(:,centerc+50:s2)=0;
    clear L num Mask centerr centerc px py mag E Bone
    %***********************************************************

    [r c]=find(V==1);
    minc=min(c);
    maxc=max(c);
    minr=min(r);
    maxr=max(r);
    clear r c 
    
    YY=Y(minr:maxr,minc:maxc);
    [ss1 ss2]=size(YY);
    
    %vertebra içinde kalan spinal kanalýn bulunmasý
    %***********************************************************
    C=ones(ss1,ss2);
    C(2:ss1-1,2:ss2-1)=0; %vertebra'yý çevreleyen pencere
    Canal=zeros(ss1,ss2);
    
    TTemp(:,:)=V(minr:maxr,minc:maxc);
    NonBone=ones(ss1,ss2);
    NonBone(TTemp==1)=0;
    
    [L num]=bwlabel(NonBone);
    for ii=1:num
        A=zeros(ss1,ss2);
        A(L==ii)=1;
        temp=A.*C;
        if ((length(nonzeros(temp))==0) && (length(nonzeros(A))>=50)) %sýnýrda olmayan ve 50 pikselden büyük olan 
                                                                    %kýsým
            mm=YY(A==1);
            mmm=mean(mm);
            clear mm

            if (mmm<145)
                Canal(A==1)=1;
            end
            clear mmm
        end 
        clear A temp
    end
    clear TTemp
    
    [L num]=bwlabel(Canal);
    
    if (num==1)
        [r c]=find(L==1);
        minc1=min(c);
        maxc1=max(c);
        minr1=min(r);
        maxr1=max(r);
        wi=maxc1-minc1+1;
        le=maxr1-minr1+1;
        if (wi/le<=2 && wi/le>=0.5)
            Canal=zeros(ss1,ss2);
            Canal(L==1)=1;
            Sp(minr:maxr,minc:maxc)=Canal(:,:);
        end
        clear r c minc1 maxc1 minr1 maxr1
    end
    clear C L num Canal NonBone minr maxr minc maxc
    %***********************************************************
    
    if (length(nonzeros(Sp))>0)
        se=[0,1,0; 1,1,1; 0,1,0];
        Sp=imerode(Sp,se);
        Sp=imfill(Sp,'holes');
        ind=i;
        Z=zeros(s1,s2);
        Z(:,:)=Z1(:,:,ind);
        Z(Sp==1)=5;
        Z1(:,:,ind)=Z(:,:);
        ftitle='results.mat';
        save(ftitle,'Z1');
        clear Z Sp inVolume Z1
        break;
    end
    clear Sp
end