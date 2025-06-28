function value=newMatrixSIR(A,eA)
[row,col]=size(A);
[erow,ecol]=size(eA);
if ((row~=erow)||(col~=ecol))
    disp('Matrix dimensions must agree!');
else
    rA=[];
    ms=0;
    tA=A; teA=eA;
    sq=sqrt(sum(tA.^2)); tA=tA./sq(ones(row,1),:);
    sq=sqrt(sum(teA.^2)); teA=teA./sq(ones(erow,1),:);
   
    for i=1:col
        [trow,tcol]=size(teA);
        wteA=[teA,-teA];
        tempai=tA(:,i);
%         for j=1:2*col
            myA=sqrt(sum((abs(tempai(:,ones(1,2*tcol))-wteA)).^2));
            [MinmyA,ind]=min(myA);
%         end
%         ms=ms+MinmyA/norm(tempai,2);
        rA=[rA,wteA(:,ind)];
        teA(:,mod(ind-1,tcol)+1)=[];  
    end
    ms=sqrt(sum((tA(:)-rA(:)).^2));
    value=-20*log10(ms/sqrt(sum((tA(:)).^2)));
%     value=-20*log10(ms/col);
end


