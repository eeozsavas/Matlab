function F=vario(I)

[s1 s2]=size(I);
k=1;

for i=1:s1
    for j=1:s2 
        F{i,j}(1)=0;
        F{i,j}(2)=0;
        F{i,j}(3)=0;
        F{i,j}(4)=0;
    end
end

%semivariogram hesaplama
for i=1+k:s1-k
    for j=1+k:s2-k

        SV0=0;
        SV45=0;
        SV90=0;
        SV135=0;
        m=0;

        W(:,:)=I(i-k:i+k,j-k:j+k); %pikselin k*k komsulugu

        for ii=1:2*k+1
            for jj=1:2*k
                SV0=SV0+(W(ii,jj)-W(ii,jj+1))^2;
                m=m+1;
            end
        end
        SV0=SV0/m;

        m=0;
        for ii=2:2*k+1
            for jj=1:2*k
                SV45=SV45+(W(ii,jj)-W(ii-1,jj+1))^2;
                m=m+1;
            end
        end
        SV45=SV45/m;

        m=0;
        for ii=2:2*k+1
            for jj=1:2*k+1
                SV90=SV90+(W(ii,jj)-W(ii-1,jj))^2;
                m=m+1;
            end
        end
        SV90=SV90/m;

        m=0;
        for ii=2:2*k+1
            for jj=2:2*k+1
                SV135=SV135+(W(ii,jj)-W(ii-1,jj-1))^2;
                m=m+1;
            end
        end
        SV135=SV135/m;

        SV(1)=SV0;
        SV(2)=SV45;
        SV(3)=SV90;
        SV(4)=SV135;

%         F{i,j}(1)=mean(SV);
%         F{i,j}(2)=std(SV);
%         F{i,j}(3)=log((SV(1)/SV(3))+(SV(3)/SV(1))+(SV(2)/SV(4))+(SV(4)/SV(2)));

        F{i,j}(1)=SV(1);
        F{i,j}(2)=SV(2);
        F{i,j}(3)=SV(3);
        F{i,j}(4)=SV(4);
        clear SV W
    end
end

%     clear I 

for m=1:4
    for i=1:s1
        for j=1:s2
            R(i,j)=F{i,j}(m);
        end
    end
    figure;imshow(R,[])
    clear R
end
%     clear F