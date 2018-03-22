%Main function-liver, kidney, spinal canal, bone segmentation

%Z1-->Background (0)
%     Patient (1)
%     Bone (2)


close all;
clear all;
clc;

se=[0,1,0; 1,1,1; 0,1,0]; %4-neighbour structuring element for morphology

%sliceNumber=125; %number of the slices
sliceNumber=1; %number of the slices

%midSlice=ceil(sliceNumber/2); %center slice
midSlice=50; %center slice

cd ('IsmailKalafat'); %directory of the slices
s1=512;
s2=512;
Z1=zeros(s1,s2,sliceNumber);
inVolume=zeros(s1,s2,sliceNumber);

%*********************************************************************
%reading the slices

%for i=1:sliceNumber
for i=50:50
    filename = sprintf('IM%01d', i); 
    X = dicomread(filename);
    X=X-1024; %rescale slope ve intercept, y=x*1+m
    inVolume(:,:,i)=X(:,:);
    
%    if (i==1) %DICOM features
     if (i==50) %DICOM features
        dInfo=dicominfo(filename);
        pSize=dInfo.PixelSpacing;
        pRegion=pSize(1)*pSize(2); %1 pixel area (mm2)
        px=round(340/pSize(1)); %patient detection için kullanilacak 340x170'lik alan
        py=round(170/pSize(2));
        px=round(px/2);
        py=round(py/2);
        sx=round(s1/2);
        sy=round(s2/2);
    end
    
    %patient delineate  
    Y=zeros(s1,s2);
    Y(X>=-175)=1; %thresholding with -175
    
    [r c]=find(Y==1);
    left=min(c);
    right=max(c);

    if (i<=midSlice/2) %for removing the table, etc...
        for jj=left:right
            temp=Y(:,jj);
            if (sum(temp)>0)
                tempp=find(temp);
                ii=tempp(length(tempp));
                sumC=sum(temp(ii-10:ii));
                if (sumC<11)
                    Y(ii-10:ii,jj)=0;
                end
            end
            clear temp tempp sumC
        end  
    end
    
    Y=imopen(Y,se); %opening three times for removing the table, etc...
    Y=imopen(Y,se); 
    Y=imopen(Y,se);
   
    overLap=zeros(s1,s2); %overlap control with 340x170 mm2 area   
    overLap(sy-py:sy+py,sx-px:sx+px)=1;

    Patient=zeros(s1,s2);
    [L num]=bwlabel(Y);
    for ii=1:num
        temp=zeros(s1,s2);
        temp(L==ii)=1;
        stemp=length(nonzeros(temp));
        temp2=overLap.*temp;
        %regions >=800 mm2 and overlap with the 340x170 mm2 area are
        %selected as patient
        if (stemp*pRegion>=800 && length(nonzeros(temp2))>0)
            Patient(L==ii)=1; %marking as patient
        end
        clear temp stemp temp2
    end
    clear L num overLap

    Patient=imclose(Patient,se);
    Patient=imclose(Patient,se);
    Patient=imclose(Patient,se);
    Patient=imfill(Patient,'holes');
    Z1(:,:,i)=Patient(:,:); 
    clear Patient r c top down left right X Y
end
cd ..; %return to the main directory


% %*********************************************************************
% % ilk kesitlerde kalan masa, vb... nesnelerin eliminasyonu için
% Y(:,:)=Z1(:,:,midSlice);
% [r c]=find(Y==1);
% down=max(r);
% clear Y r c 
% for i=1:midSlice-1    
%     Y(:,:)=Z1(:,:,i);
%     Y(down+1:s1,:)=0;
%     Z1(:,:,i)=Y(:,:);
%     clear Y
% end
% %*********************************************************************
% % son kesitlerde üst kisimda kalan masa, vb... nesnelerin eliminasyonu için
% for i=midSlice+1:sliceNumber
%     Patient=Z1(:,:,i-1);
%     [r c]=find(Patient==1);
%     top=min(r);
%     clear r c
%     
%     temp(:,:)=Z1(:,:,i);
%     [r c]=find(temp==1);
%     toptemp=min(r);
%     
%     if (toptemp<=top-5)
%         temp=temp.*Patient;
%         Z1(:,:,i)=temp(:,:);
%     end
%     clear r c temp toptemp Patient top
% end   
% fprintf('Body segmentation completed\n');


%*********************************************************************
%patient olmayan kisimlarin atilmasi
idx=find(Z1==1); 
[r c h] = ind2sub(size(Z1), idx);

minr=min(r);
maxr=max(r);
minc=min(c);
maxc=max(c);
minh=min(h); 
maxh=max(h);
Z2=Z1(minr:maxr,minc:maxc,minh:maxh); %binary volume'da patient kismi
inVolume2=inVolume(minr:maxr,minc:maxc,minh:maxh); %HU degerlerinde patient kismi
s1=size(Z2,1);
s2=size(Z2,2);
clear inVolume Z1 idx r c h


% Bone Segmentation 
% *********************************************************************
for i=1:sliceNumber
    Bone=zeros(s1,s2);
    X(:,:)=inVolume2(:,:,i);
    Body(:,:)=Z2(:,:,i);
    [r c]=find(Body~=0);
    top=min(r)+20; %patient sinirlarinin üstten ve alttan 20 piksellik kismi atiliyor
    down=max(r)-20;
    le=min(c)+10; %patient sinirlarinin sagdan ve soldan 10 piksellik kismi atiliyor
    ri=max(c)-10;
   
    X(find(Body~=1))=-10000; %background'u at
    X(1:top,:)=-10000;
    X(down:s1,:)=-10000;
    X(:,1:le)=-10000;
    X(:,ri:s2)=-10000;
    
    Bone(find(X>=145))=1;
    Bone=imfill(Bone,'holes');
    [L num]=bwlabel(Bone);
    Bone=zeros(s1,s2);
    
    for ii=1:num
        sz=length(find(L==ii));
        if ((sz*pRegion>=25))
            Bone(find(L==ii))=1;
        end
    end

    Bone=imclose(Bone,se);
    Bone=imfill(Bone,'holes');
    Body(find(Bone==1))=2;
    Z2(:,:,i)=zeros(s1,s2);
    Z2(:,:,i)=Body(:,:);
    clear X Bone Body r c top down L num tmp
end


%*********************************************************************
%convex hull of the bone structures (rib cage)
for i=1:sliceNumber
    Img(:,:)=inVolume2(:,:,i);
    
    thetas = 0:(pi/6):(pi/2); 
    fres= sqrt(2).*(2.^(0:6));
    Ga=gabor_fn(Img,thetas,fres);
    
%     F=vario(Img);
    
    Bone(:,:)=Z2(:,:,i);
    [x y]=find(Bone==2);
    k = convhull(x,y);
    BW = poly2mask(y(k), x(k), size(Bone,1),size(Bone,2)); %kemik yapilarini cevreleyen kisim-binary mask
    
    minx=min(x);
    maxx=max(x);
    miny=min(y);
    maxy=max(y);
    rightHalfImg(:,:)=Img(minx:maxx,miny:ceil((miny+maxy)/2)); %kemik yapilarini cevreleyen kismin sag yarisi-HU degerleri
    rightHalfBW(:,:)=BW(minx:maxx,miny:ceil((miny+maxy)/2)); %kemik yapilarini cevreleyen kismin sag yarisi-binary mask
    BoneHalf(:,:)=Bone(minx:maxx,miny:ceil((miny+maxy)/2)); %kemik yapilarinin 1 oldugu maskin sag yarisi-binary mask

    ImgVals=rightHalfImg(rightHalfBW==1);
    minI=min(ImgVals);
    maxI=max(ImgVals);
    rightHalfImg=round((rightHalfImg-minI)*255/(maxI-minI)); %0-255 arasina map etme
    rightHalfImg(rightHalfBW==0)=0;
    rightHalfImg(BoneHalf==2)=255;

    rightHalfImg=uint8(rightHalfImg);
    [counts,bins]=imhist(rightHalfImg);
    counts2=counts(2:254);
    bins2=bins(2:254);
    [C,ind]=max(counts2); % 0 ve 255 haric frekansi en büyük olan gray-level
    thrUp=ind+ind*1/100;
    thrDown=ind-ind*1/100;
    
    clear ImgVals minI maxI
    ImgVals=Img(BW==1);
    minI=min(ImgVals);
    maxI=max(ImgVals);
    Img=round((Img-minI)*255/(maxI-minI)); %0-255 arasina map etme
    Img(BW==0)=0;
    Img(Bone==2)=255;
    imtool(Img,[]);
    
    LiverMask=zeros(size(Img));
    LiverMask(Img>thrDown & Img<thrUp)=1;
%     figure; imshow(LiverMask,[]);
       
    LiverMask=imclose(LiverMask,se);
    LiverMask=imfill(LiverMask,'holes');
%     figure; imshow(LiverMask,[]);
    
    maxArea=1;
    comp=1;
    [L num]=bwlabel(LiverMask);     
    for ii=1:num
        cpix=length(L==ii);
        if (cpix>maxArea)
            maxArea=cpix;
            comp=ii;
        end
    end
      
    LiverMask(L~=comp)=0;
    clear L num
    figure; imshow(LiverMask,[]);
    
    fuzzyLSM(Ga,LiverMask,-0.5);
    
%     RG=fuzzyRG(Ga,LiverMask);
%     figure; imshow(RG,[]);
    
    clear x y minx maxx miny maxy Ga Bone k BW rightHalfImg rightHalfBW BoneHalf counts bins counts2 bins2 LiverMask Img RG
end

% %*********************************************************************
% %drawing contours
% for i=1:sliceNumber
%    
%     Img(:,:)=Z2(:,:,i);
%     
%     %patient contour
%     Img2=zeros(s1,s2);
%     Img2(Img~=0)=1;
%     Border = bwboundaries(Img2,8);    
%     figure; imshow(inVolume2(:,:,i),[]);
%     hold on;
%     for k = 1:numel(Border)
%          plot(Border{k}(:,2), Border{k}(:,1), 'r', 'Linewidth', 1)
%     end
%     clear Border
%     
%     %bone contour
%     Img2=zeros(s1,s2);
%     Img2(Img==2)=1;
%     Border = bwboundaries(Img2,8);    
%     for k = 1:numel(Border)
%          plot(Border{k}(:,2), Border{k}(:,1), 'y', 'Linewidth', 1)
%     end
%     clear Border Img Img2    
% end

% Bone Correction in Consecutive Slices (3D) 
% *********************************************************************
% for i=2:sliceNumber-1
%     tempZ1(:,:)=Z1(:,:,i-1);
%     L1=zeros(s1,s2);
%     L1(find(tempZ1==4))=1;
%     
%     tempZ2(:,:)=Z1(:,:,i);
%     L2=zeros(s1,s2);
%     L2(find(tempZ2==4))=1;
%     
%     tempZ3(:,:)=Z1(:,:,i+1);
%     L3=zeros(s1,s2);
%     L3(find(tempZ3==4))=1;
%     
%     for ii=1:s1
%         for jj=1:s2
%             if ((L2(ii,jj)==0) && (L1(ii,jj)==1) && (L3(ii,jj)==1))
%                 L2(ii,jj)=1;
%             end
%             if ((L2(ii,jj)==1) && (L1(ii,jj)==0) && (L3(ii,jj)==0))
%                 L2(ii,jj)=0;
%             end
%         end
%     end
%     L2=imfill(L2,'holes');
%     LT(:,:)=Z1(:,:,i);
%     LT(find(LT==4))=1;
%     LT(find(L2==1))=4;
%     Z1(:,:,i)=zeros(s1,s2);
%     Z1(:,:,i)=LT(:,:);
%     clear L1 L2 L3 tempZ1 tempZ2 tempZ3 LT
% end

% 
% %spinal kanal segmentasyonu
% ind=findFirstClosedSp(sliceNumber);
% Z1=SpExtract(ind,sliceNumber);
% fprintf('Spinal canal segmentation completed\n');
% 
% %segmentasyon sonuclarini results.mat dosyasina kaydediyor
% ftitle='results.mat';
% save(ftitle,'Z1');
% 
% clear all