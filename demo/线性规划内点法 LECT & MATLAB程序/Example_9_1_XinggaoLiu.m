function Example_9_1_XinggaoLiu
close all;
clear all;
clc;
global c
Nmax=10;
c0=16;
c=c0;
F = @def_f;
p1x=[0;5/3;2;0;0];p1y=[2.5;0;0;3;2.5];p1z=[4;4;0;0;4];
hold on
plot3(p1x,p1y,p1z);
p2x=[0;5/3;0;0];p2y=[0;0;2.5;0];p2z=[4;4;4;4];
plot3(p2x,p2y,p2z);
xlabel('x_1');ylabel('x_2');zlabel('x_3');

for i=1:7
    %x = fsolve(fun,x0,options)
    x0=[0.1;0.1;0.1];
    Y = fsolve(F,x0);
    x1(i)=Y(1);
     x2(i)=Y(2);
      x3(i)=Y(3);
    c=0.2*c;
    plot3(Y(1),Y(2),Y(3),'o')
end
x_min=Y
plot3(x1,x2,x3,'-b');
plot3(Y(1),Y(2),Y(3),'*r');
end

function F =def_f(x)
global c
F(1) = -1+3*c/(6-3*x(1)-2*x(2)-0.25*x(3))-c/x(1);
F(2) = -1+2*c/(6-3*x(1)-2*x(2)-0.25*x(3))-c/x(2);
F(3)=-5+0.25*c/(6-3*x(1)-2*x(2)-0.25*x(3))+c/(4-x(3))-c/x(3);
end



