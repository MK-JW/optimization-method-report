function Algorithm_13_7_WushengLU_Iter_0to6
% f(x) = exp(x(1)+3*x(2)-0.1) + exp(x(1)-3*x(2)-0.1) + exp(-x(1)-0.1)
f=@(x)exp(x(1)+3*x(2)-0.1) + exp(x(1)-3*x(2)-0.1) + exp(-x(1)-0.1);
g=@(x)[exp(x(1)+3*x(2)-0.1)+exp(x(1)-3*x(2)-0.1)-exp(-x(1)-0.1);3*exp(x(1)+3*x(2)-0.1)-3*exp(x(1)-3*x(2)-0.1)];
n=2;
%-  find x_star and p_star using MATLAB -------------------
[x_star,p_star] = fminsearch(f,[0.5,1]); 
%--------------------------------------------------------------------------
%--  contour near the optimal point  -------------------------------------
x0=[0.65;-0.5];A0=[1,0;0,1];
for j=1:6
subplot(2,3,j)
hold
vf=[];vf(1)=p_star;
for i=1:4
    vf(i+1) = 1.5*vf(i);
end
zoom=[-3 2.4 -1.5 1.5];
tx=zoom(1):0.05:zoom(2);ty=zoom(3):0.05:zoom(4);
[x,y]=meshgrid(tx,ty);
f = exp(x+3*y-0.1) + exp(x-3*y-0.1) + exp(-x-0.1);
contour(x,y,f,vf,'linewidth',0.8)
%---------------an ellipse at x0 -----------------------------------------

DrawEllipse(x0,A0,zoom);
%axis equal
%--------- gradient at x0 and cut plane -------------------------------
gv = g(x0);
gn = gv/norm(gv,2);
plot([x0(1),x0(1)+gn(1)],[x0(2),x0(2)+gn(2)],'-k','linewidth',1);

%-------find direction norm to grad at x0 -------
d=[gv(2);-gv(1)]./norm(gv,2);
plot([x0(1)-d(1),x0(1)+d(1)],[x0(2)-d(2),x0(2)+d(2)],'--k','linewidth',1);

%--------------an ellipse at xn -----------------------------------
gk=gv/(sqrt(gv'*A0*gv));
xn=x0-A0*gk/(n+1);
An=n^2/(n^2-1)*(A0-(2/(n+1))*A0*(gk*gk')*A0);
plot(xn(1),xn(2),'*r','linewidth',1);
DrawEllipse(xn,An,zoom);
str0=num2str(j-1);str1=num2str(j);
str=['E',str0,' & E',str1];
xlabel('x_1');ylabel('x_2');
title(str);
x0=xn;
A0=An;
plot(x_star(1),x_star(2),'*m')
end
plot(x_star(1),x_star(2),'*m')
end

function DrawEllipse(x0,A,zoom)
% zoom=[xmin xmax ymin ymax]
plot(x0(1),x0(2),'ro')
hold on
A=inv(A);
syms x y;
if nargin < 3
    h1=ezplot(A(1,1)*(x-x0(1)).^2+(A(1,2)+A(2,1))*(x-x0(1)).*(y-x0(2))+A(2,2)*(y-x0(2)).^2==1);
else
    h1=ezplot(A(1,1)*(x-x0(1)).^2+(A(1,2)+A(2,1))*(x-x0(1)).*(y-x0(2))+A(2,2)*(y-x0(2)).^2==1,zoom);
end
set(h1,'color','k','linestyle',':','linewidth',1.2); 
%axis equal
end
