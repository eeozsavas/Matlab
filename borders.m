%right-left lung separation'da cagriliyor

function S = borders(C,X)

[s1,s2]=size(C);
[r c]=find(C==1);
top=min(r);
down=max(r);
left=min(c);
right=max(c);

md=round((top+down)/2);
dst=round((right-left+1)/3);
clear r c

di=zeros(s1,s2);
di(top:md,left+dst:left+2*dst)=1-C(top:md,left+dst:left+2*dst);
[L num]=bwlabel(di);
maxArea=1;
for n=1:num
    temp=zeros(s1,s2);
    temp(L==n)=1;
    if (length(nonzeros(temp))>maxArea)
        maxArea=length(nonzeros(temp));
        ind=n;
    end
    clear temp
end

[r1 c1]=find(L==ind);
first=min(r1);
[indr,II]=min(r1);
ind=c1(II);

S=C;

k=first;
while (true)
    
    val(1)=X(k-1,ind-1);
    val(2)=X(k-1,ind);
    val(3)=X(k-1,ind+1);
    [d,i]=max(val);

    if (i==1)
        ind=ind-1;
    elseif (i==3)
        ind=ind+1;
    else
        ind=ind;
    end
    clear val
    
    if (S(k-1,ind)==0)
        break
    else
        S(k-1,ind)=0;
    end
    k=k-1;
end