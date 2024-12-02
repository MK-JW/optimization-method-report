function MainBacktrack_Armijo
close all;clear all;clc;
xk=5;
dk=1;
rho = 0.05;alpha0 = 4;
fk = fun(xk);
dphik = gfun(xk)'*dk; 
%---- plot phi(a)=f(xk+a*dk) -----------------
ae = 5;
a = 0:0.02:ae;
y = xk+a*dk;
phi=fun(y);

figure(1)
hold on
plot(a,phi,'k','LineWidth',0.9)
title('f( x_k + a*d_k) , f(x)=sin(3x)/x ,  x_k=5, d_k=1')

% LA - tanget line at xk : y = phi'(0)*a+b  s.t. fk = b  
%      pass [0,fk] and [aA,fA], fA = phi'(0)*aA + fk => aA= (fA-fk)/phi'(0)
fA =-0.2; aA = (fA-fk)/dphik;
plot([0,aA],[fk,fA],'b','LineWidth',1)
% LB - y = ruo*phi'(0)*a+b  s.t. fk = b =>  yB=ruo*phi'(0)*a+fk                                              
%      pass [0,fk] and [ae,fB]
fB=rho*dphik*ae+fk; 
plot([0,ae],[fk,fB],'g','LineWidth',1)
%-------------------------------------------------------------
[alpha]=Backtracking(@fun,@gfun,xk,dk,alpha0,rho)
%---------------------------------------------
y = xk+alpha*dk;
phi=fun(y);
plot(alpha(1),phi(1),'.r')
plot(alpha,phi,'ob')
plot(alpha(end),phi(end),'*r')
for i=1:length(alpha)
    str=num2str(i);
    text(alpha(i)+0.1,phi(i)+0.01,str)
end
str1=num2str(rho);str2=num2str(alpha0);
str=['\rho = ',str1,' ; \alpha_0 = ',str2];
text(1.5,-0.17,str)
end

function f=fun(x)
f=sin(3*x)./x;
end

function g=gfun(x)
% syms x
% f=sin(3*x)./x;
% g = diff(f,x)
g = (3*cos(3*x))./x - sin(3*x)./(x.^2);
end

function [alpha]=Backtracking(fun,gfun,xk,dk,alpha0,rho)
% alpha_star = alpha(end)  
%rho = 1.e-2;
%alpha0 = 4;
beta = 0.95;

fk = feval(fun,xk);
gk = feval(gfun,xk);

alphac = alpha0;
imax = 21;  alpha = zeros(imax,1);
alpha(1) = alpha0;
i=0;
while i < imax
    i = i + 1;
    xc = xk + alphac*dk;
    alpha(i) = alphac;
    if feval(fun,xc) < fk + rho*alphac*gk'*dk              
        alpha = alpha(1:i);
        break
    else
        alphac = beta*alphac;
    end
end
end