function example_5_1_CH05
close;
clear;
clc;
%-----------------------------------------------
%刘兴高-应用最优化方法及MATLAB实现CH05-5.1
%例5.1  目标函数f(x)=-1/(x1^2+x2^2+2);
%  x0=(2;2);  tol=1e-6
%------------------------------------------------
x_initial=[2;2];
tolerance=1e-6;
[x_optimal,f_optimal,k]=Steepest_Descent(@f_test1,@g_test1,x_initial,tolerance)
t=-2:0.1:2;
[X1,X2]=meshgrid(t);
Y=-1./(X1.^2+X2.^2+2);
contour(X1,X2,Y,'ShowText','on')
hold on
plot(x_optimal(1),x_optimal(2),'*r')
end
function f_test1=f_test1(x)
x1=x(1);
x2=x(2);
f_test1=-1/(x1^2+x2^2+2);
end
function g_test1=g_test1(x)
x1=x(1);
x2=x(2);
g1=2*x1/(x1^2+x2^2+2)^2;
g2=2*x2/(x1^2+x2^2+2)^2;
g_test1=[g1;g2];
end