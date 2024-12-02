function Example08_01_Chong
close all;clear all;clc;
x0=[4;2;-1];
kmax=21;
epsi=1.e-8;    epsi1=1.e-8;
a00=1.e-2;     a0=a00+1.e-2;
xk=x0;xr=zeros(length(x0),kmax);fr=zeros(kmax,1);
k=1;
xr = zeros(length(x0),kmax); xr(:,k)=xk; 
fr = zeros(kmax,1);          fr(k)= fun(xk);
ar = ones(kmax+1,1);         ar(1)=a00;ar(2)=a0;
gk = gfun(xk);

while norm(gk,2) > epsi && k < kmax
   dk = -gk;
   [ar(k+2),xn] = SecantMethod(@dpfun,xk,dk,a00,a0,epsi1);
   k=k+1;
   xk=xn;       gk    = gfun(xk);
   xr(:,k)=xk;  fr(k) =  fun(xk);
end
fr=fr(1:k);   ar=ar(1:k+1);   xr=xr(:,1:k);
FigureExample08_1(k,ar,xr,fr)
end

function [a,xn] = SecantMethod(fun,xk,dk,a00,a0,epsi1)
% dphi at a00,a0
kmax=10;
al=a00;ac=a0;k=1;
dphil = fun(al,xk,dk);
dphic = fun(ac,xk,dk);
while abs(dphic) > epsi1 && abs(dphic-dphil) > 1.e-3 && k < kmax
   an = ac-(ac-al)*dphic/(dphic-dphil);
   k = k+1;
   al = ac; dphil = dphic;
   ac = an; dphic = fun(ac,xk,dk);
end
a = ac; xn = xk+a*dk;
end

function f = fun(x)
f = (x(1)-4).^4+(x(2)-3).^2+4*(x(3)+5).^4;
end

function g = gfun(x)
g = [4*(x(1)-4).^3;2*(x(2)-3);16*(x(3)+5).^3];
end

function dphi = dpfun(a,xk,dk)
x = xk + a*dk;
dphi = dk'*[4*(x(1)-4).^3;2*(x(2)-3);16*(x(3)+5).^3];
end

function FigureExample08_1(k,ar,xr,fr)
%----figure step-size-----------------
figure(1)
stem(0:k,ar)
xlabel('k');ylabel('\alpha');
%---- k-x k-f ----------------
figure(2)
subplot(2,2,1)
plot(0:k-1,xr(1,1:k),'-ob')
xlabel('k');ylabel('x_1');
str1=num2str(xr(1,k));   str=['x^*_1 = ',str1]; title(str)
subplot(2,2,2)
plot(0:k-1,xr(2,1:k),'-ob')
xlabel('k');ylabel('x_2');
str2=num2str(xr(2,k));   str=['x^*_2 = ',str2];  title(str)
subplot(2,2,3)
plot(0:k-1,xr(3,1:k),'-ob')
xlabel('k');ylabel('x_3');
str3=num2str(xr(3,k));   str=['x^*_3 = ',str3];  title(str)
subplot(2,2,4)
semilogy(0:k-1,fr,'-or')
xlabel('k');ylabel('f(x)');
str0=num2str(fr(end));   str=['f^* = ',str0];   title(str)
end

