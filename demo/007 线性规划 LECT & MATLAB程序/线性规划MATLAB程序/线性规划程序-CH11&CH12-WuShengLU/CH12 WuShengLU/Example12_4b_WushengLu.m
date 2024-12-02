function Example12_4b_WushengLu
close all;
clear all;
clc;

A=[-2,2,1,-1,0;1,4,-1,0,-1];
b=[1;1];
c=[2;9;3;0;0];
tao0=5;
epdsilong=1.e-5;
[p n]=size(A);
rou=12*sqrt(n);

[w0]=example_12_4b_Constr_Var_elim;

[w_min,f_min,w,alpha,f,K]=Algorithm_12_4(A,b,c,w0,tao0,rou,epdsilong);

fileID = fopen('example12_4b_WushengLu.txt','w');
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

function [w0]=example_12_4b_Constr_Var_elim
A=[-2,2,1,-1,0;1,4,-1,0,-1];
b=[1;1];
c=[2;9;3;0;0];
[p n]=size(A);
[U,S,V]= svd(A);
%AA=U*S*V'
Vr=V(:,(p+1):n);
b_p=pinv(A)*b;
Vpha=[1,0;0,1;0,1];
VR=Vr*Vpha;
pha1=0.8;%Let pha=[0.8;0.5];
pha0=0.5;
pha=[pha1;pha0];
x0=VR*pha+b_p;
lamda0=[1;1];
miu0=c-A'*lamda0;
w0=[x0;lamda0;miu0];
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

Nmax=50;

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
%delta_lamda=Y*A*y;
delta_lamda=linsolve(A*D*A',A*y);
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
