function example_5_2_CH05
close;
clear;
clc;
%-----------------------------------------------
%刘兴高-应用最优化方法及MATLAB实现CH05-5.2
%例5.2  目标函数f(x)=x1^2+x2^2+x1*x2+2;
%  x0=(1;-4);  tol=1e-6
%------------------------------------------------
x_initial=[1;-4];
tolerance=1e-6;
[x_optimal,f_optimal,k]=Steepest_Descent(@f_test2,@g_test2,x_initial,tolerance)
t=-5:0.1:5;
[X1,X2]=meshgrid(t);
Y=X1.^2+X2.^2+X1.*X2+2;
contour(X1,X2,Y,'ShowText','on')
hold on
plot(x_initial(1),x_initial(2),'or')
plot(x_optimal(1),x_optimal(2),'*r')
end
function f_test2=f_test2(x)
x1=x(1);
x2=x(2);
f_test2=x1^2+x2^2+x1*x2+2;
end
function g_test2=g_test2(x)
x1=x(1);
x2=x(2);
g1=2*x1+x2;
g2=2*x2+x1;
g_test2=[g1;g2];
end