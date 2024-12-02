function Example_13_7_WushengLU
close all;clear all; clc;
% Algorithm 13.8 Ellipsoid method for CP constrained problems
% Example13.7
% f = (x1-5*x2+4).^2+(7*x1+11*x2-18).^4
% c1 =  2*x1-x2- 6 >= 0
% c2 = -2*x1-x2+10 >= 0
% c3 = -3*x1+x2+15 >= 0
% c4 =  3*x1+x2- 9 >= 0
 x0=[4;-1];
 A0=[36,0;0,9];
 maxk=5000;
 epsilon=1e-7;
 k=1; 
%
format long
%---------- find p_star using MATLAB ---------------------------------
%[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
% A*x <= b
A=[-2,2,3,-3;1,1,-1,-1]';b=[-6;10;15;-9];
%f=(x1-5*x2+4).^2+(7*x1+11*x2-18).^4;
[x_star,p_star] = fmincon(@fun,x0,A,b,[],[],[],[],[]);
str1=num2str(x_star(1));str2=num2str(x_star(2));str3=num2str(p_star);
str=['MATLAB sovled ', ' ; x* = [ ',str1,'; ',str2,' ] ; p* = ',str3];
disp(str);
%--  end p_star -------------------------------------------------------
n=length(x0);
xc=x0;Ac=A0;
%---------------------------------------------------------------------
fc=fun(xc);
str1=num2str(xc(1));str2=num2str(xc(2));str3=num2str(fc);
str=['Ellipsoid method Begin...... ', ' ; x0 = [ ',str1,'; ',str2,' ] ; f0 = ',str3];
disp(str);
%---------------------------------------------------------------------

%-----  begin Figure (1) for initial case  ---------------------------
figure(1)
hold on 
zoom=[2.5 5.5 -3.5 2.5];
% plot the constrants
h1 = ezplot(@(x1,x2) 2*x1-x2-6,zoom);
set(h1,'color','k','linestyle','--','linewidth',0.9);
h2 = ezplot(@(x1,x2) -2*x1-x2+10,zoom);
set(h2,'color','k','linestyle','--','linewidth',0.9);
h3 = ezplot(@(x1,x2) -3*x1+x2+15,zoom);
set(h3,'color','k','linestyle','--','linewidth',0.9);
h4 = ezplot(@(x1,x2) 3*x1+x2- 9,zoom);
set(h4,'color','k','linestyle','--','linewidth',0.9);
text(3.5,1.5,'c_1');text(3.5,-1.9,'c_4');
text(4.5,1.5,'c_2');text(4.5,-1.9,'c_3');
title('min f(x) = (x1-5*x2+4).^2+(7*x1+11*x2-18).^4');
plot(xc(1),xc(2),'or');% initial point
% conotur f
fc=fun(xc);
dv=(fc-70)/10;
vf=70:dv:fc;
[x1,x2]=meshgrid(2.5:0.05:5.5,-3.5:0.05:2.5);
f=(x1-5*x2+4).^2+(7*x1+11*x2-18).^4;
contour(x1,x2,f,vf)
% plot the initial ellipse
DrawEllipse(xc,Ac,zoom);
xlabel('x_1');ylabel('x_2');
%--------end figure(1) for the initial case -----------------------------

% record iterations
xr=zeros(n,maxk);fr=zeros(maxk,1);
xr(:,1)=xc;fr(1)=fc;
%----------------------------------------------------------
while k < maxk
    cval = cfun(xc);
    [cmin,Ic]=min(cval);
    if cmin < 0
        hk = hfun(xc,Ic);
        gk1=hk/sqrt(hk'*Ac*hk);
        if cmin+sqrt(hk'*Ac*hk) < 0
           break;
        end
        xn = xc-Ac*gk1/(n+1);
        An = n^2/(n^2-1)*(Ac-2/(n+1)*Ac*(gk1*gk1')*Ac);

    else
        gk     = gfun(xc);
        gk1    = gk/sqrt(gk'*Ac*gk);
        gammak = sqrt(gk'*Ac*gk);
        if(gammak < epsilon)        
         break;
        end
        xn = xc-Ac*gk1/(n+1);
        An = n^2/(n^2-1)*(Ac-2/(n+1)*Ac*(gk1*gk1')*Ac);
    end
    k=k+1;
    xc=xn;
    Ac=An;
    %--- record record iterations
    xr(:,k)=xc;fr(k)=fun(xc);
    %-------------------------------
end

%--- record record iterations
    xr=xr(:,1:k);fr=fr(1:k);
%-------------------------------
%-------begin figure(2) for iterations and final ellipse -----------------
figure(2)
hold on 
zoom=[1.6 4.5 -1.5 0.7];
% plot the constrants
h1 = ezplot(@(x1,x2) 2*x1-x2-6,zoom);
set(h1,'color','k','linestyle','-','linewidth',0.8);
h4 = ezplot(@(x1,x2) 3*x1+x2- 9,zoom);
set(h4,'color','k','linestyle','-','linewidth',0.8);
text(3.2,0.3,'c_1');text(3.2,-1.2,'c_4');

% iteration points
plot(xr(1,1),xr(2,1),'or');% initial point
plot(xr(1,:),xr(2,:),'--k','linewidth',0.8);
plot(xr(1,:),xr(2,:),'.b','MarkerSize',10);
plot(xr(1,end),xr(2,end),'*r','linewidth',0.9);
% conotur f
zoom=[1.6 4.5 -1.5 0.7];
vf=fr(1:min(6,length(fr)));
vf=[vf;fr(end)];
[x1,x2]=meshgrid(1.6:0.05:4.5,-1.5:0.05:0.7);
f=(x1-5*x2+4).^2+(7*x1+11*x2-18).^4;
contour(x1,x2,f,vf,'linewidth',0.9)
% plot the final ellipse
DrawEllipse(xc,Ac,zoom);
xlabel('x_1');ylabel('x_2');
title('min f(x) = (x1-5*x2+4).^2+(7*x1+11*x2-18).^4');
%--------  end  figure(2) ------------------------------

%--------- figure (3) iter v.s. f(x_k) - p* -----------------
figure(3)
semilogy(1:k,abs(fr-p_star),'b','linewidth',0.9);
xlabel('k');ylabel('| f(x^{(k)}) - p^* |');
%-------------------------------------------------------------

fc=fun(xc);
str0=num2str(k);str1=num2str(xc(1));str2=num2str(xc(2));str3=num2str(fc);
str=['Ellipsoid method Completed.  Iter = ', str0,' ; x* = [ ',str1,'; ',str2,' ] ; f* = ',str3];
disp(str);
end

function DrawEllipse(x0,A,zoom)
plot(x0(1),x0(2),'r+')
hold on
A=inv(A);
syms x y;
if nargin < 3
    h1=ezplot(A(1,1)*(x-x0(1)).^2+(A(1,2)+A(2,1))*(x-x0(1)).*(y-x0(2))+A(2,2)*(y-x0(2)).^2==1);
else
    h1=ezplot(A(1,1)*(x-x0(1)).^2+(A(1,2)+A(2,1))*(x-x0(1)).*(y-x0(2))+A(2,2)*(y-x0(2)).^2==1,zoom);
end
set(h1,'color','b','linestyle',':','linewidth',0.9); 
%axis equal
end



function fval = fun(x)
fval = (x(1)-5*x(2)+4).^2+(7*x(1)+11*x(2)-18).^4;
end

function gval = gfun(x)
gval=[2*(x(1)-5*x(2)+4)+28*(7*x(1)+11*x(2)-18)^3;
    -10*(x(1)-5*x(2)+4)+44*(7*x(1)+11*x(2)-18)^3];
end

function cval = cfun(x)
cval = [ 2*x(1)-x(2)- 6;
        -2*x(1)-x(2)+10;
        -3*x(1)+x(2)+15;
         3*x(1)+x(2)- 9];
end

function hval = hfun(x,i)
% c1 =  2*x1-x2- 6 >= 0
% c2 = -2*x1-x2+10 >= 0
% c3 = -3*x1+x2+15 >= 0
% c4 =  3*x1+x2- 9 >= 0
% h =[-▽c1,-▽c2,-▽c3,-▽c4]
hval=-[  2, -2, -3,  3;
        -1, -1,  1,  1];
if nargin < 2
    %hval = hval;
else
    hval = hval(:,i);
end
end


