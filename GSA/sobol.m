function [S,ST]=sobol(D,nPop,VarMin,VarMax,myFunction)
M=D*2;%

VarMin=[VarMin,VarMin];
VarMax=[VarMax,VarMax];
p= sobolset(M);
% https://www.cnblogs.com/zhubinglong/p/12260292.html

r=p(1:nPop,:);
R=VarMin+r.*(VarMax-VarMin);

A=R(:,1:D);%
B=R(:,D+1:end);
AB=zeros(nPop,D,D);
for i=1:D
    tempA=A;
    tempA(:,i)=B(:,i);
    AB(1:nPop,1:D,i)=tempA;
end

YA=zeros(nPop,1);%
YB=zeros(nPop,1);
YAB=zeros(nPop,D);%
YA=myFunction(A(:,:));
YB=myFunction(B(:,:));
for j=1:D
    YAB(:,j)=myFunction(AB(:,:,j));
end

%myFunction(A(1:5,:))
%find(~isreal(YAB(:,3)))

VarX=zeros(D,1);%
S=zeros(D,1);

VarY=var([YA;YB],1,'omitnan');%
for i=1:D
    for j=1:nPop
        VarX(i)=VarX(i)+YB(j).*(YAB(j,i)-YA(j));
    end
    VarX(i)=1./nPop.*VarX(i);%
    S(i)=VarX(i)./VarY;
end

%%
EX=zeros(D,1);
ST=zeros(D,1);
for i=1:D
    for j=1:nPop
        EX(i)=EX(i)+(YA(j)-YAB(j,i))^2;
    end
    EX(i)=1/(2*nPop)* EX(i);%
    ST(i)=EX(i)/VarY;
end
end