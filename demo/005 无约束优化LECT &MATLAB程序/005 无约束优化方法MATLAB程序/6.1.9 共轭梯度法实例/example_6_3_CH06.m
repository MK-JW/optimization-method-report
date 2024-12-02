function example_6_3_CH06
close;
clear;
clc;
%-----------------------------------------------
%刘兴高-应用最优化方法及MATLAB实现CH06-6.3
%例6.1  目标函数f(x)=1+x1+x2+x3+x4+x1*x2+x1*x3+x1*x4+x2*x3+x2*x4+x3*x4+x1^2+x2^2+x3^2+x4^2-0.4*exp(-x5^2-6*x6^2);
%  x0=(-4;0;-4;-1;1;1);  tol=1e-6
%------------------------------------------------
x_initial=[-4;0;-4;-1;1;1];
tolerance=1e-6;
[x_optimal,f_optimal,k]=Conjugate_gradient_DY(@f_test3,@g_test3,x_initial,tolerance)
end
function f_test3=f_test3(x)
x1=x(1);x2=x(2);x3=x(3);x4=x(4);x5=x(5);x6=x(6);
f_test3=1+x1+x2+x3+x4+x1*x2+x1*x3+x1*x4+x2*x3+x2*x4+x3*x4+x1^2+x2^2+x3^2+x4^2-0.4*exp(-x5^2-6*x6^2);
end
function g_test3=g_test3(x)
x1=x(1);x2=x(2);x3=x(3);x4=x(4);x5=x(5);x6=x(6);
g1=2*x1+x2+x3+x4+1;g2=x1+2*x2+x3+x4+1;
g3=x1+x2+2*x3+x4+1;g4=x1+x2+x3+2*x4+1;
g5=(4*x5)/(5*exp(x5^2+6*x6^2));g6=(24*x6)/(5*exp(x5^2+6*x6^2));
g_test3=[g1;g2;g3;g4;g5;g6];
end