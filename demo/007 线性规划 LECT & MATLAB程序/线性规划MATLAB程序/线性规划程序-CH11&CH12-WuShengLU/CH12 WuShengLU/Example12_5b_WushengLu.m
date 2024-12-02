function Example12_5b_WushengLu
close all;
clear all;
clc;

A=[-2,2,1,-1,0;1,4,-1,0,-1];
b=[1;1];
c=[2;9;3;0;0];

epsi=1.e-8;
[p n]=size(A);
rou=12*sqrt(n);

%Start_point='Auto';
Start_point='Give';
if Start_point=='Auto'
    w0=[];%if we want to compute a starting point,let 'w0=[] '
else
    x0=[1;0.1;0.1;2;5];
    lmd0=[-1;1];
    mu0=[1.0;0.1;0.2;1;10];
    w0=[x0;lmd0;mu0];
end

[w_min,f_min,w,alpha,f,K]=Algorithm_12_5(A,b,c,w0,rou,epsi);

fileID = fopen('example12_5b_WushengLu.txt','w');
fprintf(fileID,'%9s %32s %30s %18s %34s\n','alpha','x','f','lamda','miu');

for k=1:K
   fprintf(fileID,'%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\r\n',alpha(k),w(1:5,k),f(k),w(6:12,k));
end
fclose(fileID);

%=============plot k-x1 k-x2 k-x3 k-x4 k-x5 k-f=========
subplot(2,3,1)
kv=1:K;
plot(kv,w(1,:))
title('Subplot 1: x_1')

subplot(2,3,2)
plot(kv,w(2,:))
title('Subplot 2: x_2')

subplot(2,3,3)
plot(kv,w(3,:))
title('Subplot 3: x_3')

subplot(2,3,4)
plot(kv,w(4,:))
title('Subplot 4: x_4')

subplot(2,3,5)
plot(kv,w(5,:))
title('Subplot 5: x_5')

subplot(2,3,6)
plot(kv,f)
title('Subplot 6: f')
%============================
end

function [w0]=Initial_point(A,b,c)
[p n]=size(A);
e=ones(n,1);
nl0=n+1;
nl1=n+p;
nm0=n+p+1;
nm1=n+n+p;
AAi=inv(A*A');
x_b=A'*AAi*b;
lam_b=AAi*A*c;
mu_b=c-A'*lam_b;
d_x=max(-1.5*(min(x_b)),0);
d_mu=max(-1.5*(min(mu_b)),0);
x_h=x_b+d_x*e;
mu_h=mu_b+d_mu*e;
d_x_h=0.5*((x_h'*mu_h)/(e'*mu_h));
d_mu_h=0.5*((x_h'*mu_h)/(e'*x_h));
w0=[x_h+d_x_h*e;lam_b;mu_h+d_mu_h*e];
pp=A*w0(1:n)-b;%if pp=0, w0(1:n) is feasible point
end

function [w_min,f_min,w,alpha,f,K]=Algorithm_12_5(A,b,c,w0,rou,epsi)
%=====initialization==================
if isempty(w0)
    w0=Initial_point(A,b,c)
end
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
%tao=tao0;

Nmax=50;

while duality_gap>epsi & k<Nmax
    k=k+1;
    [w_k,alpha_k,f_k]=Algorithm_12_5_w(A,b,c,w0,rou);
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

function [w_k,alpha_k,f_k]=Algorithm_12_5_w(A,b,c,w0,rou)
%===============================================================
%===============================================================
[p n]=size(A);
x=w0(1:n);
lamda=w0(n+1:n+p);
miu=w0(n+p+1:2*n+p);

e=ones(n,1);

tao=(miu'*x)/(n+rou);

M_inv=diag(1./miu);
D=diag(x./miu);  %   D=M_inv*X;  X=diag(x);
Y=inv(A*D*A');
%y=x-tao*M_inv*e;
y=x-tao./miu;

r_p=b-A*x;
r_d=c-A'*lamda-miu;

delta_lamda=Y*(A*y+A*D*r_d+r_p);
%delta_lamda=linsolve(Y,(A*y+A*D*r_d+r_p));
delta_miu=-A'*delta_lamda+r_d;
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