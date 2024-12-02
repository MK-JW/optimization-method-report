function example12_3_WushengLu
close all;
clear all;
clc;
format long
A=[1 1 1];b=1;c=[-2;1;-3];
x00=1/3;x0=[x00;x00;x00];
n=length(x0);
tao_0=0.1;sigma=0.1;tao=tao_0;
epsilong_inn=0.0001;epsilong_out=1.e-4;
%Step 1
l=1;
x(:,1)=x0;
norm_x=1;norm_ad=1;
f(1)=c'*x0;
a_L1=0;%for GoldenSection alpha_low
while norm_x>=epsilong_out
    tao_1=1/tao;
    k=1;
    y0=x(:,l);
    y00=y0;
    norm_ad=1;   
    while norm_ad>=epsilong_inn & k<50
       %========search direction=======
        X=diag(y0);
        X2=diag(y0.*y0);
        AA=A*X2*A';
        bb=A*(X2*c-tao*y0);
        lamda=AA\bb;
        d=y0+tao_1*X2*(A'*lamda-c);d=d/norm(d);
        %================================
        %=======alpha_hat=========
        jj=0;
        for j=1:n
           if d(j)<0
              jj=jj+1;
              if jj==1
                  a_min=-y0(j)/d(j);
              else
                 a_min=min(a_min,-y0(j)/d(j));
              end
           end
        end
        a_U1=0.99*a_min;
        %for GoldenSection alpha_upper
       %======================================================
       %=======Golden-section for alpha=======================
       [as,fs] = golden_sect_alpha(y0,d,tao,a_L1,a_U1,1e-10);
       %=======================================================         
         y=y0+as*d;
         y0=y;
         k=k+1;
         a=as;
         norm_ad=norm(a*d);
    end
    l=l+1;
    x(:,l)=y;
    norm_x=norm(y-y00);
    f(l)=fs;
    tao=sigma*tao;
end
format long
x
f
x_MIN=x(:,l);
%============constrant plane================
xc=[1;0;0;1];yc=[0;1;0;0];zc=[0;0;1;0];
plot3(xc,yc,zc);
%====================================================
%============plot central path==================
hold on
plot3(x(1,:),x(2,:),x(3,:),'-om')
[mm,nn]=size(x);
plot3(x(1,nn),x(2,nn),x(3,nn),'*r')
f_min=min(f);
str=num2str(f_min);
st=strcat('f_min=  ',str);
title(st)
xlabel('x1');
ylabel('x2');
zlabel('x3');
%============================================
end

function fv=fun_def(x,d,a,tao)
c=[-2;1;-3];
fv=c'*(x+a*d)-tao*(log(x(1)+a*d(1))+log(x(2)+a*d(2))+log(x(3)+a*d(3)));
end

function [am,fm] = golden_sect_alpha(x,d,tao,aL,aU,rou)
if aU>aL
    r=0.5*(3-sqrt(5));
    a0=aL;b0=aU;
    Nmax=500;k=0;
    L=b0-a0;
    while L>=rou & k<Nmax
        k=k+1;
        LG=r*L;
        a1=a0+LG;
        b1=b0-LG;
        fa1=fun_def(x,d,a1,tao);
        fb1=fun_def(x,d,b1,tao);
        if fa1<fb1
            b0=b1;
        else
            a0=a1;
        end
        L=b0-a0;
    end
    am=0.5*(a1+b1);
    fm=fun_def(x,d,am,tao);
end
end