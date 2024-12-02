function Example_13_6_Iter_1TO4
%Example 13.6 (1) 
close;clear;clc;
gname='gex';fname='fex';
epsi=1e-7;
x0=[4,-1]';
A0=[36,0;0,9];
n=length(x0);
x=x0;A=A0;

xr=[];yr=[];
xr=[xr x(1)];yr=[yr x(2)];
AA0=A0;xx0=x0;
%==================Iter 0=========================================
N=12;
for kr=1:N
%subplot(2,2,1)
figure(kr)
[x1,x2]=meshgrid(-8:0.05:12);
fx = (x1-5*x2+4).^2+(7*x1+11*x2-18).^4;
v=1:300:1000;
contour(x1,x2,fx,v,'k')
hold on
h1=ezplot(@(xt,yt)ellipse(xt,yt,xx0,AA0),[-8,12,-8,12]);
set(h1,'color','b','linestyle',':'); % E_k-1
%
h2=ezplot(@(xt,yt)ellipse(xt,yt,x,A),[-8,12,-8,12]);
set(h2,'color','b','linestyle','--'); %E_k
%-------------------------------------------------------------
xx0=x;
AA0=A;
gk = gex(x);
gak = sqrt(gk'*A*gk);
gkt = gk/gak;
plot(xr(1),yr(1),'or')
plot(xr(1:kr),yr(1:kr),'c')
plot(xr(1:kr),yr(1:kr),'.r')
plot(xr(kr),yr(kr),'*m')
hold on
beta=2;
gk1=x(1)+gk(1)/(sqrt(gk(1)^2+gk(2)^2));
gk2=x(2)+gk(2)/(sqrt(gk(1)^2+gk(2)^2));
plot([x(1),beta*gk1],[x(2),beta*gk2],'b')
hold on
beta=2;
gkt1=x(1)+gkt(1)/(sqrt(gkt(1)^2+gkt(2)^2));
gkt2=x(2)+gkt(2)/(sqrt(gkt(1)^2+gkt(2)^2));
plot([x(1),beta*gkt1],[x(2),beta*gkt2],'-.g')
hold on
p = ezplot(@(xt,yt)hyperplane(xt,yt,x,gk),[-7,11,-7,11]);
set(p,'color','b','linestyle','-.'); % hyperplane
nums=num2str(kr-1);
str=['Iteration  ' nums];
title(str);
axis equal

%-------   x_next  and A_next  --------------------------------------
x = x - A*gkt/(n+1);
xr=[xr x(1)];yr=[yr x(2)];
Az = A*gkt;      
A = n^2*(A - 2*(Az*Az')/(n+1))/(n^2-1);
%------------------------------------------------------------
end
end



function z = hyperplane(x,y,x0,g)
z=g(1)*(x-x0(1))+g(2)*(y-x0(2));
end

%=============   ellipse  ===================================
function z = ellipse(x,y,x0,A)
aa=A(1,1);bb=A(1,2);cc=A(2,1);dd=A(2,2);
t=aa*dd-cc*bb;
a=dd/t;d=aa/t;b=-bb/t;c=-cc/t;
x1=x0(1);y1=x0(2);
z=a*(x-x1).^2+(b+c)*(x-x1).*(y-y1)+d*(y-y1).^2-1;
end

% ============= Evaluate a subgradient at x =============
% f(x) = (x1 − 5*x2 + 4)^2 + (7*x1 + 11*x2 − 18)^4
 function g = gex(x)
 x1 = x(1)-5*x(2)+4; 
 x2= (7*x(1) + 11*x(2)-18)^3;
 g1=  2*x1 + 28*x2;
 g2=-10*x2 + 44*x2;
 g=[g1;g2];
 end
% ========================================================

% ========================================================
% f(x) = (x1−5*x2+4)^2+(7*x1+11*x2−18)^4
function f = fex(x)
f = (x(1)-5*x(2)+4)^2 + (7*x(1)+11*x(2)-18)^4;
end
% =====================
% 