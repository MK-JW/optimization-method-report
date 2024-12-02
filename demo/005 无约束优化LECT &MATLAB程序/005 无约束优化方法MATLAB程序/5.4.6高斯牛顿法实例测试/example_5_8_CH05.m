function example_5_8_CH05
close;
clear;
clc;
%-----------------------------------------------
%刘兴高-应用最优化方法及MATLAB实现CH05-5.8
%例5.8  目标函数F(x)=(x1+5)^2+(x2+8)^2+(x3+7)^2+2*x1^2*x2^2+4*x1^2*x3^2;
%  [f1;f2;f3;f4;f5]=[x1+5;x2+8;x3+7;2^0.5*x1*x2;2*x1*x3]
%  x0=(4;13;1);  tol=1e-6
%------------------------------------------------
x_initial=[4;13;1];
tolerance=1e-6;
[x_optimal,F_optimal,k]=Guass_Newton(@F_test2,@G_test2,@J_test2,x_initial,tolerance)
end

function F_test2=F_test2(x)
x1=x(1);
x2=x(2);
x3=x(3);
F_test2=(x1+5)^2+(x2+8)^2+(x3+7)^2+2*(x1^2)*(x2^2)+4*(x1^2)*(x3^2);
end
function G_test2=G_test2(x)
x1=x(1);
x2=x(2);
x3=x(3);
G1=4*x1*(x2^2)+8*x1*(x3^2)+2*x1+10;
G2=4*x2*(x1^2)+2*x2+16;
G3=8*x3*(x1^2)+2*x3+14;
G_test2=[G1;G2;G3];
end
function J_test2=J_test2(x)
x1=x(1);
x2=x(2);
x3=x(3);
J11=1;
J12=0;
J13=0;
J21=0;
J22=1;
J23=0;
J31=0;
J32=0;
J33=1;
J41=2^(1/2)*x2;
J42=2^(1/2)*x1;
J43=0;
J51=2*x3;
J52=0;
J53=2*x1;
J_test2=[J11,J12,J13;J21,J22,J23;J31,J32,J33;J41,J42,J43;J51,J52,J53];
end