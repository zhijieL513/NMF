function [eA]=kverticefinding(V,A0)
[m,n]=size(A0);
T=length(V);
enA=A0;
myp=[];
myset=combnk(1:n,m-1);
[myrow,mycol]=size(myset);
[irow,icol]=size(myset);
for i=1:irow
    tmpeA=enA(:,myset(i,:));
    [tmpv,tmpd]=eig(tmpeA*tmpeA');
    myp(:,i)=tmpv(:,1);
end

oldmyp=myp;

for i=1:myrow %K
    for t=1:10
        nc=mean(enA,2); normc=sqrt(sum(nc.^2));  nc=nc/normc;
        mysign=nc'*myp(:,i);
        tmpind=find(mysign*myp(:,i)'*V<=0);
        T0=round(.05*T/n);
        if  (T0<=m-1)
            T0=m;
        end
        if (length(tmpind)>=T0)
            tmpM=V(:,tmpind);
            [tmpv,tmpd]=eig(tmpM*tmpM');
            myp(:,i)=tmpv(:,1);
        end
    end
end


for i=1:myrow
    corr(i)=abs(myp(:,i)'*oldmyp(:,i));
    if corr(i)<0.97
        myp(:,i)=oldmyp(:,i);
    end
end

enA=[];
for i=1:myrow
    mysign=zeros(myrow,1);
    for j=1:mycol
        mysign=mysign|(myset(:,j)==i);
    end
    tmpind=find(mysign==1);
    tmpM=myp(:,tmpind);
    [tmpv,tmpd]=eig(tmpM*tmpM');
    enA(:,i)=(tmpv(:,1));
end
eA=enA;