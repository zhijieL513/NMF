%
% Code for paper "Conic Hull Fitting Based Dicionary Matrix Learning for Nonnegative Matrix Factorization"
% Author: Zhijie Lin.
% 

m=4;
n=3;
K=n;
T=1000; 

ourSNR=[]; GVavasisSNR=[]; combinedSIR=[];
ourVOL=[]; GVavasisVOL=[]; 

for seed=1:100
    rand('state',seed); A=rand(m,n);
    sqa=sqrt(sum(A.^2)); A=A./sqa(ones(m,1),:);
    rand('state',seed+1); S=rand(n,T);

    X=A*S;

    SNR=20;
    rand('state',555); 
    noise=rand(m,T);
    sn=mean(noise.^2,2);
    sX=mean(X.^2,2); 
    myc=sX./(sn*(10^(SNR/10)));
    X=X+diag(sqrt(myc))*noise;

    [tU,tD]=eig(X*X'); ttX=tD(m-n+1:m,m-n+1:m)^(-.5)*tU(:,m-n+1:m)'*X;

    myset=combnk(1:n,n-1);
    [myrow,mycol]=size(myset);

    [J,normM,U] = FastSepNMF(X,K,1); enA=X(:,J);
    
    [combinedeA]=kverticefinding(ttX,tD(m-n+1:m,m-n+1:m)^(-.5)*tU(:,m-n+1:m)'*enA);
    sqea=sqrt(sum(combinedeA.^2)); combinedeA=combinedeA./sqea(ones(n,1),:);
    combinedSIR(seed)=newMatrixSIR(A,tU(:,m-n+1:m)*tD(m-n+1:m,m-n+1:m)^(.5)*combinedeA);
end % seed

mean(combinedSIR)