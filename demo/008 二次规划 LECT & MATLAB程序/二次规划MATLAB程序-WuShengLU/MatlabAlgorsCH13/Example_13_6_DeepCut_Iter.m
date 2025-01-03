function Example_13_6_DeepCut_Iter
%Example 13.6 (1) WushengLU
close;clear;clc;

DEEPCUT=1; % DEEPCUT=1 deep cut. DEEPCUT=0 Classical Ellisoid method
epsi=1e-7;
%x0=[4,-1]';
%A0=[36,0;0,9];
x0=[31,11]';A0=[1600,0;0,256];
n=length(x0);
x=x0;A=A0;

xr=[];yr=[];f_save = []; f_best = []; 
xr=[xr x(1)];yr=[yr x(2)];
AA0=A0;xx0=x0;
%==================Iter=========================================
N=6;
for kr=1:N
    %----------------- contour --------------------
%subplot(2,3,kr)
figure(kr)
[x1,x2]=meshgrid(-8:0.05:12);
fx = (x1-5*x2+4).^2+(7*x1+11*x2-18).^4;
v=1:300:1000;
contour(x1,x2,fx,v,'k')
%-----------------ellipsoid at x_last----------------------------
hold on
h1=ezplot(@(xt,yt)ellipse(xt,yt,xx0,AA0),[-8,12,-8,12]);
set(h1,'color','b','linestyle',':'); % E_k-1
%-----------------ellipsoid at x_current----------------------------
h2=ezplot(@(xt,yt)ellipse(xt,yt,x,A),[-8,12,-8,12]);
set(h2,'color','b','linestyle','--'); %E_k
%-------------------------------------------------------------
xx0=x;
AA0=A;
    f = fex(x);% find function value at current x
    g = gex(x);% subgradient at current x
	% convergence test
	if sqrt(g'*A*g) <= epsi
        xs=x;fs=f;
		disp('converged')
		break
	end
    % save current and best function values 
    f_save = [f_save f];
    f_best = [f_best min(f_save)];
    h = f - f_best(end);
gk = g;
gak = sqrt(gk'*A*gk);
gkt = gk/gak;
plot(xr(1),yr(1),'or')
plot(xr(1:kr),yr(1:kr),'c')
plot(xr(1:kr),yr(1:kr),'.r')
plot(xr(kr),yr(kr),'*m')
%---------------plot grad ------------------------
   hold on
   beta=2;% plot grad at x
   gk1=x(1)+gk(1)/(sqrt(gk(1)^2+gk(2)^2));
   gk2=x(2)+gk(2)/(sqrt(gk(1)^2+gk(2)^2));
   plot([x(1),beta*gk1],[x(2),beta*gk2],'b')
   hold on
   beta=2;% plot A*grad at x
   gkt1=x(1)+gkt(1)/(sqrt(gkt(1)^2+gkt(2)^2));
   gkt2=x(2)+gkt(2)/(sqrt(gkt(1)^2+gkt(2)^2));
   plot([x(1),beta*gkt1],[x(2),beta*gkt2],'-.g')
%-------------plot the hyperplane-----------------
   hold on
   p = ezplot(@(xt,yt)hyperplane(xt,yt,x,gk,0),[-7,11,-7,11]);
   set(p,'color','b','linestyle',':'); % hyperplane
   if DEEPCUT > 0 && h > 0
       p = ezplot(@(xt,yt)hyperplane(xt,yt,x,gk,h),[-7,11,-7,11]);
       set(p,'color','b','linestyle','-.'); % hyperplane
   end
   
   nums=num2str(kr-1);
   str=['Iteration  ' nums];
   title(str);
   axis equal
%-------------update--x_next  and A_next------------------------------
   h = f - f_best(end);
    if DEEPCUT > 0 && h > 0
		% update with deep cut.
		alpha = h/sqrt(gk'*A*gk);
		gkt = gk/sqrt(gk'*A*gk);
		x = x - (1+alpha*n)/(n+1)*A*gkt;
		A = n^2*(1-alpha^2)/(n^2-1)*(A-2*(1+alpha*n)/((n+1)*(1+alpha))*A*(gkt*gkt')*A);
        xr=[xr x(1)];yr=[yr x(2)];
	else
		% update ellipsoid with shallow cut.
		gkt = g/sqrt(gk'*A*gk);
		x = x - 1/(n+1)*A*gkt;
		A = n^2/(n^2-1)*(A-2/(n+1)*A*(gkt*gkt')*A);
        xr=[xr x(1)];yr=[yr x(2)];
    end 
    print -deps Example_13_6_DeepCut_Iter.eps
end
end



function z = hyperplane(x,y,x0,g,h)
z=g(1)*(x-x0(1))+g(2)*(y-x0(2))+h;
end

%=============   ellipse  ===================================
function z = ellipse(x,y,x0,A)
%aa=A(1,1);bb=A(1,2);cc=A(2,1);dd=A(2,2);
t=A(1,1)*A(2,2)-A(2,1)*A(1,2);
a=A(2,2)/t;d=A(1,1)/t;b=-A(1,2)/t;c=-A(2,1)/t;
%x1=x0(1);y1=x0(2);
z=a*(x-x0(1)).^2+(b+c)*(x-x0(1)).*(y-x0(2))+d*(y-x0(2)).^2-1;
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