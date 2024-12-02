%Draw an ellipse E={x| (x-x0)'*inv(A)*(x-x0) = 1}
function DrawAnEllipse
x0=[4;-1];
A0=[36,0;0,9];
zoom=[2.5 5.5 -3.5 2.5];
DrawEllipse(x0,A0,zoom)
hold on
DrawEllipse(x0,A0)
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