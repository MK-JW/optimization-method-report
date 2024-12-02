function Example12_4a_WushengLu
close all;
clear all;
clc;
A=[1 1 1];
b=1;
c=[-2 1 -3]';
tao0=5;
epdsilong=1.e-6;
[p n]=size(A);
rou=7*sqrt(n);

x0=[0.34450640;0.28549374;0.36999986];
lamda0=-16.51351863;
miu0=[14.51351863;17.51351863;13.51351863];
w0=[x0;lamda0;miu0];

[w_min,f_min,w,alpha,f,K]=Algorithm_12_4(A,b,c,w0,tao0,rou,epdsilong);


fileID = fopen('example12_4a_WushengLu.txt','w');
fprintf(fileID,'%9s %28s %18s %18s %18s\n','alpha','x','f','lamda','miu');

for k=1:K
 fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',alpha(k),w(1,k),w(2,k),w(3,k),f(k),w(4,k),w(5,k),w(6,k),w(7,k));
end
fclose(fileID);
%============constrant plane================
xc=[1;0;0;1];yc=[0;1;0;0];zc=[0;0;1;0];
plot3(xc,yc,zc);
%====================================================
%============plot central path==================
hold on
plot3(w(1,:),w(2,:),w(3,:),'-om')
plot3(w(1,K),w(2,K),w(3,K),'*r')

str=num2str(f_min);
st=strcat('f_min=  ',str);
title(st)
xlabel('x1');
ylabel('x2');
zlabel('x3');
%============================================
end
function [w_min,f_min,w,alpha,f,K]=Algorithm_12_4(A,b,c,w0,tao0,rou,epdsilong)
[p n]=size(A);
nl0=n+1;
nl1=n+p;
nm0=n+p+1;
nm1=n+n+p;
%---------------------------------------------------
%-                        w(1:n,    k)=x(k)        -
%-  nl0=n+1;nl1=n+p-1;    w(nl0:nl1,k)=lamda(k)    -
%-  nm0=n+p+1;nm1=n+n+p;  w(nm0:nm1,k)=miu(k)      -
%---------------------------------------------------
k=1;
duality_gap=w0(nm0:nm1)'*w0(1:n);
f(k)=c'*w0(1:n);
alpha(k)=0;
w(:,k)=w0;
tao=tao0;

Nmax=1000;

while duality_gap>epdsilong & k<Nmax
    k=k+1;
    [w_k,alpha_k,f_k]=Algorithm_12_4_w(A,c,w0,rou);
    w(:,k)=w_k;
    alpha(k)=alpha_k;
    f(k)=f_k;
    duality_gap=w_k(nm0:nm1)'*w_k(1:n);
    w0=w_k;    
end
format short
K=k
format long
w_min=w(:,k)
f_min=f(k)
end

function [w_k,alpha_k,f_k]=Algorithm_12_4_w(A,c,w0,rou)
%===============================================================
%===============================================================
[p n]=size(A);
x=w0(1:n);
miu=w0(n+p+1:2*n+p);
e=ones(n,1);
tao=(miu'*x)/(n+rou);
M_inv=diag(1./miu);
X=diag(x);
%D=M_inv*X;
D=diag(x./miu);
Y=inv(A*D*A');
y=x-tao*M_inv*e;
delta_lamda=Y*A*y;
delta_miu=-A'*delta_lamda;
delta_x=-y-D*delta_miu;

a_p=1;
a_d=1;
for j=1:n
   if delta_x(j)<0
       a_p_temp=-x(j)./delta_x(j);
       if a_p>a_p_temp
           a_p=a_p_temp;
       end
   end
   if delta_miu(j)<0
       a_d_temp=-miu(j)./delta_miu(j);
       if a_d>a_d_temp
           a_d=a_d_temp;
       end
   end
end
alpha_k=0.999999*min(a_p,a_d);
delta_w=[delta_x;delta_lamda;delta_miu];
w_k=w0+alpha_k*delta_w;
f_k=c'*w_k(1:n);
end
