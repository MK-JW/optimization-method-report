function Example12_6a_WushengLu
close all;
clear all;
clc;
A=[1 1 1];
b=1;
c=[-2 1 -3]';
[p n]=size(A);

epsi=1.e-6;

%Start_point='Auto';
Start_point='Give';
if Start_point=='Auto'
    w0=[];%if we want to compute a starting point,let 'w0=[] '
else
    x0  =[0.4;0.3;0.4];
    lam0=0.5;
    mu0 =[1.0;0.5;1.0];
    w0  =[x0;lam0;mu0];
end

[w_min,f_min,w,f,K]=Algorithm_12_6(A,b,c,w0,epsi);

fileID = fopen('example12_6a_WushengLu.txt','w');
fprintf(fileID,'%22s %18s %12s %18s\n','x','f','lamda','miu');
for k=1:K
   fprintf(fileID,'%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\r\n',w(1:3,k),f(k),w(4:7,k));
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

function [w0]=Initial_point(A,b,c)
%A=[1 1 1];
%b=1;
%c=[-2 1 -3]';

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

function [w_min,f_min,w,f,K]=Algorithm_12_6(A,b,c,w0,epsi)

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
x=w0(1:n);
mu=w0(nm0:nm1);
lam=w0(nl0:nl1);
gap=mu'*x;
f(k)=c'*x;
w(:,k)=w0;
e=ones(n,1);
Nmax=50;

while gap>epsi & k<Nmax
    k=k+1;
    tau_h=gap/n;
    %===  Step 3 predictor direction 
    mui=1./mu;
    D = diag(x./mu);
    Y = inv(A*D*A');
    r_d = c - A'*lam - mu;
    d_aff_lam=Y*(b+A*D*r_d);
    d_aff_mu =r_d-A'*d_aff_lam;
    d_aff_x  =-x-D*d_aff_mu;
%=== Step4 compute tao_aff and tao_k_+_1
    a_p=1;
    a_d=1;
    for j=1:n
        if d_aff_x(j)<0
            a_p_temp=-x(j)./d_aff_x(j);
            if a_p>a_p_temp
                a_p=a_p_temp;
            end
        end
        if d_aff_mu(j)<0
            a_d_temp=-mu(j)./d_aff_mu(j);
            if a_d>a_d_temp
               a_d=a_d_temp;
            end
        end
    end
    a_aff_p = min(a_p);
    a_aff_d = min(a_d);
    
    tau_aff=[(mu+a_aff_d*d_aff_mu)'*(x+a_aff_p*d_aff_x)]/n;
    sigma_k=(tau_aff/tau_h)^3;
    tau=sigma_k*tau_h;

    %========Step 5 Compute corrector direction delta_c_w============
    d_X=diag(1./d_aff_x);
    y=mui.*d_aff_x.*d_aff_mu-tau*mui;
    d_c_lam=Y*A*y;
    d_c_mu =-A'*d_c_lam;
    d_c_x  =-y-D*d_c_mu;
%=================================================================

%======Step 6 Obtain search direction=============================
   d_lam=d_c_lam+d_aff_lam;
   d_mu =d_c_mu  +d_aff_mu;
   d_x  =d_c_x    +d_aff_x;
%===  a_kp and a_kd=====
   a_kp=1/0.99;
   a_kd=1/0.99;
   for j=1:n
       if d_x(j)<0
          a_kp_temp=-x(j)/d_x(j);
          if a_kp>a_kp_temp
              a_kp=a_kp_temp;
          end
       end
       if d_mu(j)<0
          a_kd_temp=-mu(j)/d_mu(j);
          if a_kd>a_kd_temp
              a_kd=a_kd_temp;
          end
       end
   end
   a_kp=min(0.99*a_kp,1);
   a_kd=min(0.99*a_kd,1);
%======Step 7 =================
   x_k=x+a_kp*d_x;
   lam_k=lam+a_kd*d_lam;
   mu_k=mu+a_kd*d_mu;
   f_k=c'*x_k;  
   w(:,k)=[x_k;lam_k;mu_k];
   f(k)=f_k;
   gap=mu_k'*x_k;
   
   x=x_k;
   lam=lam_k;
   mu=mu_k;
end
format short
K=k
format long
w_min=w(:,k)
f_min=f(k)
end

