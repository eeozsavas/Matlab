% spinal kanal segmentasyonu

function Z1 =SpExtract(ind,sliceNumber)

se=[0,1,0; 1,1,1; 0,1,0]; %4-neighbour structuring element for morphology

slc_ctrl=zeros(sliceNumber,1);
Z1=zeros(512,512,sliceNumber);
inVolume=zeros(512,512,sliceNumber);

ftitle='HU.mat';
load(ftitle);

ftitle='results.mat';
load(ftitle);

for i=ind-1:-1:1
    slc_ctrl(i)=0;
    
    XX(:,:)=inVolume(:,:,i);

    Z(:,:)=Z1(:,:,i);
    
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
    
    clear L num Mask centerr centerc px py mag E Bone minc maxc minr maxr
    %***********************************************************

    [r c]=find(V==1);
    minc=min(c);
    maxc=max(c);
    minr=min(r);
    maxr=max(r);
    clear r c 
    
    clear YY
    YY(:,:)=Y(minr:maxr,minc:maxc);
    VRoi=zeros(s1,s2);
    VRoi(minr:maxr,minc:maxc)=1;
    [ss1 ss2]=size(YY);
    
    %vertebra içinde kalan spinal kanalýn bulunmasý
    %***********************************************************
    C=ones(ss1,ss2);
    C(2:ss1-1,2:ss2-1)=0; %vertebra'yý çevreleyen pencere
    CC=zeros(s1,s2);
    CC(minr:maxr,minc:maxc)=C(:,:);
    Canal=zeros(ss1,ss2);
    
    TTemp(:,:)=V(minr:maxr,minc:maxc);
    NonBone=ones(ss1,ss2);
    NonBone(YY>=145)=0;
    
    [L num]=bwlabel(NonBone);
    for ii=1:num
        A=zeros(ss1,ss2);
        A(L==ii)=1;
        temp=A.*C;
        if ((length(nonzeros(temp))==0) && (length(nonzeros(A))>=50)) %sýnýrda olmayan ve 50 pikselden büyük olan 
                                                                    %kýsým
            mm=(YY(A==1));
            mmm=mean(mm);
            clear mm

            if (mmm<145)
                Canal(A==1)=1;
            end
            clear mmm
        end
        clear A temp
    end
    clear TTemp L num temp
    
    [L num]=bwlabel(Canal);
    SpBefore1(:,:)=Z1(:,:,i+1);
    SpBefore=zeros(s1,s2);
    finCanal=zeros(s1,s2);
    SpBefore(SpBefore1==5)=1;
    clear SpBefore1
    for ii=1:num
        CanalS=zeros(s1,s2);
        Canal=zeros(ss1,ss2);
        Canal(L==ii)=1;
        CanalS(minr:maxr,minc:maxc)=Canal(:,:);
        ctemp=CanalS.*SpBefore;
        temp=length(nonzeros(ctemp));
        if (temp>(length(nonzeros(Canal))/2))
            finCanal(CanalS==1)=1;
        end
        clear temp CanalS Canal ctemp
    end
 
    comp=1;
    maxArea=0;
    clear L num
    [L num]=bwlabel(finCanal);
    for ii=1:num
        Canal=zeros(s1,s2);
        Canal(L==ii)=1;
        if (length(nonzeros(Canal))>maxArea)
            maxArea=length(nonzeros(Canal));
            comp=ii;
        end
        clear Canal
    end
      
    Sp=zeros(s1,s2); 
    
    if (maxArea>0)
        Sp(L==comp)=1;
        Sp=imfill(Sp,'holes');
        
        [r c]=find(Sp==1);
        minc1=min(c);
        maxc1=max(c);
        minr1=min(r);
        maxr1=max(r);
        clear r c 
        
        [r c]=find(SpBefore==1);
        minc2=min(c);
        maxc2=max(c);
        minr2=min(r);
        maxr2=max(r);
        clear r c
        
        if ((abs(minc1-minc2)>5) || (abs(maxc1-maxc2)>5) || (abs(minr1-minr2)>5) || (abs(maxr1-maxr2)>5))
            Sp=zeros(s1,s2);
        end
    end
   
    Z(:,:)=Z1(:,:,i);
    Z(Sp==1)=5;
    Z1(:,:,i)=Z(:,:);
    clear Z L num Canal Canal1 NonBone minr maxr minc maxc minr1 maxr1 minc1 maxc1 minr2 maxr2 minc2 maxc2 SpBefore
    %***********************************************************
    
    if (length(nonzeros(Sp))==0)
        
        slc_ctrl(i)=1;
        
        SpBefore1(:,:)=Z1(:,:,i+1);
        SpBefore(:,:)=zeros(s1,s2); 
        SpBefore(SpBefore1==5)=1;
        clear SpBefore1
        
        YBefore(:,:)=inVolume(:,:,i+1);
                
        [r c]=find(SpBefore==1);        

        minc=min(c);
        maxc=max(c);
        minr=min(r);
        maxr=max(r);
        centerr=round((minr+maxr)/2);
        centerc=round((minc+maxc)/2);
        clear r c

        Sptemp=zeros(s1,s2);
        vvv=23;
        svv=15;
        m=0;
       
        for i1=centerr-10:centerr+10
            for j1=centerc-10:centerc+10
                if (Y(i1,j1)>=vvv-3*svv && Y(i1,j1)<=vvv+3*svv)
                    Sptemp(i1,j1)=1;
                    m=m+1;
                end
            end
        end

        if (m==0)
            Sptemp(centerr,centerc)=1;
        end
        
        mea=23;
        st=15;

        
        FM=exp(double((-1.*((Y-mea).^2)./((st^2)*2))));
        FM(VRoi==0)=0; 
        FM(Sptemp==1)=1;

        marker=zeros(s1,s2);
        marker(Sptemp==1)=1.0;        
       
        res=imreconstruct(marker,FM); 
        thr=0.4;
        res(res<thr)=0;
        res(res>=thr)=1;
        res=imfill(res,'holes');
        
        maxArea=1;
        comp=1;
        [L num]=bwlabel(res);     
        for ii=1:num
            A=zeros(s1,s2);
            A(L==ii)=1;
            rtemp=A.*SpBefore;
            temp=length(nonzeros(rtemp));
            if (temp>maxArea)
                maxArea=temp;
                comp=ii;
            end
            clear A rtemp temp
        end
        res=zeros(s1,s2); 
        res(L==comp)=1;
        clear L num
        
        J=zeros(s1,s2);
        A=res.*CC;
        if (length(nonzeros(A))==0)
            J=res;
        else
            A=res.*SpBefore;
            J(A==1)=1;
        end
        clear A res
        J=imfill(J,'holes');
        clear Sp temp Border r c e d z a b alpha
        Z(:,:)=Z1(:,:,i);
        Z(J==1)=5;
        Z1(:,:,i)=Z(:,:);
        clear Z       
    end
    clear C CC Y YY Sp SpBefore YBefore finalvals mea st J marker Sptemp Bone VRoi vv vvv svv V FM
end


for i=ind+1:sliceNumber
    slc_ctrl(i)=0;
    XX(:,:)=inVolume(:,:,i);
    Z(:,:)=Z1(:,:,i);
    
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
    
    clear L num Mask centerr centerc px py mag E Bone minc maxc minr maxr
    %***********************************************************

    [r c]=find(V==1);
    minc=min(c);
    maxc=max(c);
    minr=min(r);
    maxr=max(r);
    clear r c 
    
    clear YY
    YY(:,:)=Y(minr:maxr,minc:maxc);
    VRoi=zeros(s1,s2);
    VRoi(minr:maxr,minc:maxc)=1;
    [ss1 ss2]=size(YY);
    
    %vertebra içinde kalan spinal kanalýn bulunmasý
    %***********************************************************
    C=ones(ss1,ss2);
    C(2:ss1-1,2:ss2-1)=0; %vertebra'yý çevreleyen pencere
    CC=zeros(s1,s2);
    CC(minr:maxr,minc:maxc)=C(:,:);
    Canal=zeros(ss1,ss2);
    
    TTemp(:,:)=V(minr:maxr,minc:maxc);
    NonBone=ones(ss1,ss2);
    NonBone(YY>=145)=0;
    
    [L num]=bwlabel(NonBone);
    for ii=1:num
        A=zeros(ss1,ss2);
        A(L==ii)=1;
        temp=A.*C;
        if ((length(nonzeros(temp))==0) && (length(nonzeros(A))>=50)) %sýnýrda olmayan ve 50 pikselden büyük olan 
                                                                    %kýsým
            mm=(YY(A==1));
            mmm=mean(mm);
            clear mm

            if (mmm<145)
                Canal(A==1)=1;
            end
            clear mmm
        end
        clear A temp
    end
    clear TTemp L num temp
    
    [L num]=bwlabel(Canal);
    SpBefore1(:,:)=Z1(:,:,i-1);
    SpBefore=zeros(s1,s2);
    finCanal=zeros(s1,s2);
    SpBefore(SpBefore1==5)=1;
    clear SpBefore1
    for ii=1:num
        CanalS=zeros(s1,s2);
        Canal=zeros(ss1,ss2);
        Canal(L==ii)=1;
        CanalS(minr:maxr,minc:maxc)=Canal(:,:);
        ctemp=CanalS.*SpBefore;
        temp=length(nonzeros(ctemp));
        if (temp>(length(nonzeros(Canal))/2))
            finCanal(CanalS==1)=1;
        end
        clear temp CanalS Canal ctemp
    end
 
    comp=1;
    maxArea=0;
    clear L num
    [L num]=bwlabel(finCanal);
    for ii=1:num
        Canal=zeros(s1,s2);
        Canal(L==ii)=1;
        if (length(nonzeros(Canal))>maxArea)
            maxArea=length(nonzeros(Canal));
            comp=ii;
        end
        clear Canal
    end
      
    Sp=zeros(s1,s2); 
    
    if (maxArea>0)
        Sp(L==comp)=1;
        Sp=imfill(Sp,'holes');
        
        [r c]=find(Sp==1);
        minc1=min(c);
        maxc1=max(c);
        minr1=min(r);
        maxr1=max(r);
        clear r c 
        
        [r c]=find(SpBefore==1);
        minc2=min(c);
        maxc2=max(c);
        minr2=min(r);
        maxr2=max(r);
        clear r c
        
        if ((abs(minc1-minc2)>5) || (abs(maxc1-maxc2)>5) || (abs(minr1-minr2)>5) || (abs(maxr1-maxr2)>5))
            Sp=zeros(s1,s2);
        end
    end
   
    Z(:,:)=Z1(:,:,i);
    Z(Sp==1)=5;
    Z1(:,:,i)=Z(:,:);
    clear Z L num Canal Canal1 NonBone minr maxr minc maxc minr1 maxr1 minc1 maxc1 minr2 maxr2 minc2 maxc2 SpBefore
    %***********************************************************
    
    if (length(nonzeros(Sp))==0)
        
        slc_ctrl(i)=1;
        
        SpBefore1(:,:)=Z1(:,:,i-1);
        SpBefore(:,:)=zeros(s1,s2); 
        SpBefore(SpBefore1==5)=1;
        clear SpBefore1
     
        YBefore(:,:)=inVolume(:,:,i-1);
                
        [r c]=find(SpBefore==1);        

        minc=min(c);
        maxc=max(c);
        minr=min(r);
        maxr=max(r);
        centerr=round((minr+maxr)/2);
        centerc=round((minc+maxc)/2);
        clear r c

        Sptemp=zeros(s1,s2);

        vvv=23;
        svv=15;
        m=0;
       
        for i1=centerr-10:centerr+10
            for j1=centerc-10:centerc+10
                if (Y(i1,j1)>=vvv-3*svv && Y(i1,j1)<=vvv+3*svv)
                    Sptemp(i1,j1)=1;
                    m=m+1;              
                end
            end
        end

        if (m==0)
            Sptemp(centerr,centerc)=1;
        end
        
        mea=23;
        st=15;
        
        FM=exp(double((-1.*((Y-mea).^2)./((st^2)*2))));
        FM(VRoi==0)=0; 
        FM(Sptemp==1)=1;

        marker=zeros(s1,s2);
        marker(Sptemp==1)=1.0;        
       
        res=imreconstruct(marker,FM); 
        thr=0.4;
        res(res<thr)=0;
        res(res>=thr)=1;
        res=imfill(res,'holes');
        
        maxArea=1;
        comp=1;
        [L num]=bwlabel(res);     
        for ii=1:num
            A=zeros(s1,s2);
            A(L==ii)=1;
            rtemp=A.*SpBefore;
            temp=length(nonzeros(rtemp));
            if (temp>maxArea)
                maxArea=temp;
                comp=ii;
            end
            clear A rtemp temp
        end
        res=zeros(s1,s2); 
        res(L==comp)=1;
        clear L num
        
        J=zeros(s1,s2);
        A=res.*CC;
        if (length(nonzeros(A))==0)
            J=res;
        else
            A=res.*SpBefore;
            J(A==1)=1;
        end
        clear A res
        J=imfill(J,'holes');
        
        clear Sp temp Border r c e d z a b alpha
        Z(:,:)=Z1(:,:,i);
        Z(J==1)=5;
        Z1(:,:,i)=Z(:,:);
        clear Z       
    end
    clear C CC Y YY Sp SpBefore YBefore finalvals mea st J marker Sptemp Bone VRoi vv vvv svv V FM    
end

for i=2:sliceNumber-1
    S11=zeros(s1,s2);
    S22=zeros(s1,s2);
    S33=zeros(s1,s2);
    S11(:,:)=Z1(:,:,i-1);
    S22(:,:)=Z1(:,:,i);
    S33(:,:)=Z1(:,:,i+1);
    S1=zeros(s1,s2);
    S2=zeros(s1,s2);
    S3=zeros(s1,s2);
    S1(S11==5)=1;
    S2(S22==5)=1;
    S3(S33==5)=1;
    for ii=1:s1
        for jj=1:s2
            if ((S2(ii,jj)==0) && (S1(ii,jj)==1) && (S3(ii,jj)==1))
                S2(ii,jj)=1;
            end
            if ((S2(ii,jj)==1) && (S1(ii,jj)==0) && (S3(ii,jj)==0))
                S2(ii,jj)=0;
            end
        end
    end
    S2=imopen(S2,se);
    S2=imfill(S2,'holes');
    maxArea=1;
    comp=1;
    [L num]=bwlabel(S2);     
    for ii=1:num
        temp=length(find(L==ii));
        if (temp>maxArea)
            maxArea=temp;
            comp=ii;
        end
        clear temp
    end
    S2=zeros(s1,s2);
    S2(L==comp)=1;
    clear L num
    S22(S22==5)=1;
    S22(S2==1)=5;
    Z1(:,:,i)=zeros(s1,s2);
    Z1(:,:,i)=S22(:,:);
    
    clear S1 S2 S3 S11 S22 S33
end

for i=1:sliceNumber
    Sp=zeros(s1,s2);
    temp(:,:)=Z1(:,:,i);
   
    X(:,:)=inVolume(:,:,i);

    Sp(temp==5)=1;
    Sp=imopen(Sp,se);
    
    if (slc_ctrl(i)==1)
  
        Border = edge(Sp);
        [r c]=find(Border==1);
        e(:,1)=r(:);
        e(:,2)=c(:);

        [z, a, b, alpha,errctr] = fitellipse(e, 'linear');
        if (errctr==0)
            d=plotellipse(z, a, b, alpha);

            d=round(d);

            Sp=zeros(s1,s2);

            for ii=1:200
                Sp(d(1,ii),d(2,ii))=1;
            end

            Sp=imfill(Sp,'holes');
            Sp(X>=145)=0;
            Sp=imfill(Sp,'holes');
            Sp=imopen(Sp,se);
        end
        
    else
        Sp=imerode(Sp,se);
    end
    
    maxArea=1;
    comp=1;
    [L num]=bwlabel(Sp);     
    for ii=1:num
        temp1=length(find(L==ii));
        if (temp1>maxArea)
            maxArea=temp1;
            comp=ii;
        end
        clear temp1
    end
    Sp=zeros(s1,s2);
    Sp(L==comp)=1;
    clear L num comp
    
    temp(temp==5)=1;
    temp(Sp==1)=5;
    Z1(:,:,i)=temp(:,:);
    clear Sp temp Border X r c e d z a b alpha
end
clear inVolume