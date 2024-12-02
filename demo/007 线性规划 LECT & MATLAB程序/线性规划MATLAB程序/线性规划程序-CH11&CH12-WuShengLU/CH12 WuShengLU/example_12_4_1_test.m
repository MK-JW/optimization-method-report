function example_12_4_1_test
close all;
clear all;
clc;
x0=[];
lamda0=0.0;
miu0=[];
w0=[x0;lamda0;miu0];
tao0=5;
epdsilong=1.e-6;
n=length(x0);
A=[1 1 1];
b=1;
c=[-2 -3 2]';
%--------------------------------------------
miu=[1;2;4];
x=[0.2;0.2;0.2];
A=[1 1 1];
[p n]=size(A);
n=3;
rou=12.5;

[w_k,alpha_k,f_k]=NPL_w(A,c,w_k_1,rou)
%--------------------------------------------
%=============================
%===   w(nx0:nx1,:)=x;     ===
%===   w(nl0:nl1,:)=lamda; ===
%===   w(nm0:nm1,:)=miu;   ===
nx0=1;nx1=n;
nlamda=length(lamda);
nl0=n+1;nl1=nl0+nlamda-1;
nm0=nl1+1;nm1=nm0+n-1;
%=============================

end

function [w,alpha,f]=NPL_w(A,c,w0,rou)
%miu=[1;2;4];
%lamda=1;
%x=[0.2;0.2;0.2];
%w0=[x;lamda;miu];
%A=[1 1 1];
%rou=12.5;
%c=[-3;-2;1]
[p n]=size(A);
x=w0(1:n);
miu=w0(n+p+1:2*n+p);
e=ones(n,1);
tao=(miu'*x)/(n+rou);
M_inv=diag(1./miu);
X=diag(x);
D=M_inv*X
Y=inv(A*D*A')
y=x-tao*M_inv*e
delta_lamda=Y*A*y
delta_miu=-A'*delta_lamda
delta_x=-y-D*delta_miu
a_p=1000;a_d=1000;
for i=1:n
    if delta_x(i)<0
        a_p=min(a_p,-x(i)/delta_x(i));
    end
    if delta_miu(i)<0
        a_d=min(a_d,-miu(i)/delta_miu(i));
    end
end
a_min=min(a_d,a_p);
alpha=0.999999*a_min;
w=w0+alpha*w0;
f=c'*w(1:n);
end
