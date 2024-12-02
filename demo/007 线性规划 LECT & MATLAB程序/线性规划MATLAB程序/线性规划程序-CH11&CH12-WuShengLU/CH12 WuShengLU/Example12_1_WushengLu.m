function Example12_1_WushengLu
close all;
clear all;
clc;
a=1;b=1;c=1;d=1;C=[-2 1 -3]';
A=[1 1 1];
%syms lamda
%eqn=a*lamda^3+b*lamda^2+c*lamda+d==0;
tao_0=5;
tao_f=0.0001;
N=20;
dtao=(tao_f-tao_0)/N;
tao=tao_0;
k=0;
fileID = fopen('example12_1_WushengLu.txt','w');
fprintf(fileID,'%9s %28s %18s %18s %18s\n','tao','x','f','lamda','miu');
for i=1:N+1
    a=1/tao;
    b=4*a+3;
    c=a+8;
    d=1-6*a;
    x=[a;b;c;d];
    s_lamda=shengjin(x);
    lamda=s_lamda(s_lamda<-3);
    %if lamda+3<1.e-8
       % test=-3-lamda
        k=k+1;
        x1(k)=tao/(-2-lamda);
        x2(k)=tao/(1-lamda);
        x3(k)=tao/(-3-lamda);
       miu=C-A'*lamda;
       f(k)=-2*x1(k)+x2(k)-3*x3(k);
    %end
    t(i)=tao;
    %fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f\n',t(i),x1(k),x2(k),x3(k),f(k));
    fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',t(i),x1(k),x2(k),x3(k),f(k),lamda,miu(1),miu(2),miu(3));
    tao=tao+dtao;
end
fclose(fileID);
%============constrant surface================
x_max=max(x1);y_max=max(x2);z_max=max(x3);
x_min=min(x1);y_min=min(x2);z_min=min(x3);
Mx=20;My=20;
dx=(x_max-x_min)/Mx;dy=(y_max-y_min)/My;
tx=x_min:dx:x_max;ty=y_min:dy:y_max;
[X Y]=meshgrid(tx,ty);
Z=1-X-Y;
mesh(X,Y,Z);
%====================================================
%============plot central path==================
hold on
plot3(x1,x2,x3,'-om')
plot3(x1(N+1),x2(k),x3(k),'*r')
f_min=min(f);
str=num2str(f_min);
st=strcat('f_m_i_n=  ',str);
title(st)
xlabel('x1');
ylabel('x2');
zlabel('x3');
%============================================
end

function realroot=shengjin(x)
%ShengJin's formula for solving cubic equations
         A=x(2).^2-3*x(1).*x(3);
         B=x(2).*x(3)-9*x(1).*x(4);
         C=x(3).^2-3*x(2).*x(4);
         theta=B.^2-4*A.*C;
         if A==0&B==0
            xroot=-x(3)./x(2);
         elseif theta>0
            Y1=A.*x(2)+3*x(1).*(-B+theta.^(1/2))/2;
            Y2=A.*x(2)+3*x(1).*(-B-theta.^(1/2))/2;
            xroot=(-x(2)-(nthroot(Y1,3)+nthroot(Y2,3)))./(3*x(1));
         elseif theta==0
            xroot(1)=-x(2)./x(1)+B./A;
            xroot(2)=-B./(2*A);
         elseif theta<0
             T=(2*A.*x(2)-3*x(1).*B)./(2*A.^(3/2));
             phi=acos(T);
             xroot(1)=-(x(2)+2*A.^(1/2).*cos(phi/3))./(3*x(1));
             xroot(2)=(-x(2)+A.^(1/2).*(cos(phi/3)+3^(1/2)*sin(phi/3)))./(3*x(1));
             xroot(3)=(-x(2)+A.^(1/2).*(cos(phi/3)-3^(1/2)*sin(phi/3)))./(3*x(1));
         end
         realroot=xroot;%only real toot
end