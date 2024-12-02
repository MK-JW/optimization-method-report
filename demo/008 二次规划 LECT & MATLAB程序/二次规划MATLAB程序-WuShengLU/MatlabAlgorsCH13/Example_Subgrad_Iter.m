function Example_Subgrad_Iter
%------------------------------------------
%   f(x) = |x|+2|y|; f(x) = |x|+4|y|
%-------------------------------------------
close;clear;clc;
gname='gex';fname='fex';
epsi=1e-7;
x0=[4,-1]';
A0=[16,0;0,16];
n=length(x0);
x=x0;A=A0;

xr=[];yr=[];
xr=[xr x(1)];yr=[yr x(2)];
AA0=A0;xx0=x0;
%==================Iter 0=========================================
N=6;
for kr=1:N
subplot(2,3,kr)
%figure(kr)
[x1,x2]=meshgrid(-8:0.05:9);
%fx=abs(x1)+2*abs(x2);
fx=abs(x1)+4*abs(x2);
v=1:4:30;
contour(x1,x2,fx,v,':k')
hold on
h1=ezplot(@(xt,yt)ellipse(xt,yt,xx0,AA0),[-8,9,-8,9]);
set(h1,'color','b','linestyle','--'); % E_k-1
%
h2=ezplot(@(xt,yt)ellipse(xt,yt,x,A),[-8,9,-8,9]);
set(h2,'color','b','linestyle','-'); %E_k
%-------------------------------------------------------------
xx0=x;
AA0=A;
gk = gex(x);
gak = sqrt(gk'*A*gk);
gkt = gk/gak;
plot(xr(1),yr(1),'or')
plot(xr(1:kr),yr(1:kr),'c')
plot(xr(1:kr),yr(1:kr),'.r')
plot(xr(kr),yr(kr),'.m')
plot(0,0,'.m')
hold on
beta=2;
gk1=x(1)+gk(1)/(sqrt(gk(1)^2+gk(2)^2));
gk2=x(2)+gk(2)/(sqrt(gk(1)^2+gk(2)^2));
plot([x(1),beta*gk1],[x(2),beta*gk2],'k')
hold on
beta=2;
gkt1=x(1)+gkt(1)/(sqrt(gkt(1)^2+gkt(2)^2));
gkt2=x(2)+gkt(2)/(sqrt(gkt(1)^2+gkt(2)^2));
plot([x(1),beta*gkt1],[x(2),beta*gkt2],':g')
hold on
p = ezplot(@(xt,yt)hyperplane(xt,yt,x,gk),[-7,9,-7,9]);
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


% ========================================================
% f(x) = |x|+2|y| ; f(x) = |x|+4|y|
function f = fex(x)
%f = abs(x(1)) + 2*abs(x(2));
f = abs(x(1)) + 4*abs(x(2));
end
% =====================
% 

function g = gex(x0)
%---------------------------------------------
%   f(x)=|x|+2|y|
%---------------------------------------------
n=length(x0);
g=zeros(n,1);
x=x0(1);
y=x0(2);
tx=1.1;% -2 < tx < 2
ty=0.8;% -1 < ty <1
if x == 0
    if y == 0
        g(1)=1;% any direction
        g(2)=-1;
    else
        g(1)=tx;
        g(2)=sign(y);
    end
else
    if y== 0
        g(1)=2*sign(x);
        g(2)=ty;
    else
        g(1)=sign(x);
        g(2)=sign(y);
    end
end
end
