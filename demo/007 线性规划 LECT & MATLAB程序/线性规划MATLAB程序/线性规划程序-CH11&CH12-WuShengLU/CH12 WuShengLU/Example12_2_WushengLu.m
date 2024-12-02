close all;
clear all;
clc;
format long
x0=[0.2;0.7;1;1;1];
N=length(x0);
A=[-2,2,1,-1,0;1,4,-1,0,-1];
b=[1;1];
c=[2;9;3;0;0];
gamma=0.9999;
epsilong=0.000001;
f0=c'*x0;
e=epsilong+1;
k=0;
fileID = fopen('example12_2_WushengLu.txt','w');
fprintf(fileID,'%9s %28s %18s\n','k','x','f');

while e>epsilong
    k=k+1;
    X=diag(x0);
    X2=diag(x0.*x0);
    AX2=A*X2;
    d=-[X2-AX2'*(inv(AX2*A'))*AX2]*c;
    i0=0;
    for i=1:N
        if d(i)<0
            i0=i0+1;
            t(i0)=-x0(i)/d(i);
        end
    end
    alpha=gamma*min(t);
    x1=x0+alpha*d;
    f1=c'*x1;
    e=abs(f1-f0)/max(1,abs(f0));
    fprintf(fileID,'%4d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',k,x1',f1,e);
    f0=f1;
    x0=x1;
end