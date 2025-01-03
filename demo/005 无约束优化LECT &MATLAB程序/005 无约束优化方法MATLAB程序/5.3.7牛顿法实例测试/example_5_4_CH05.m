function example_5_4_CH05
close;
clear;
clc;
%-----------------------------------------------
%刘兴高-应用最优化方法及MATLAB实现CH05-5.1
%例5.4  目标函数f(x)=-1/(x1^2+x2^2+2);
%  x0=(2;2);  tol=1e-6
%------------------------------------------------
x_initial=[2;2];
tolerance=1e-6;
[x_optimal,f_optimal,k]=Newton(@f_test1,@g_test1,@H_test1,x_initial,tolerance)
t=-3:0.1:3;
[X1,X2]=meshgrid(t);
Y=-1./(X1.^2+X2.^2+2);
%---------------------------------------------------------------------
%读取记录迭代点的文件，并画出迭代点
%---------------------------------------------------------------------
ex=importdata('testdata.txt');
data=ex.data;
[m,n]=size(data);
K=data(1:m,1)+1;%迭代次数
X=data(1:m,2:n-1);%迭代点
F=data(1:m,n);%迭代点的函数值
%-----------画等高线（值为F的）---------------
contour(X1,X2,Y,F)
%-------------------------------------------
%----------画迭代点轨迹---------------------
hold on
plot(X(1:m,1),X(1:m,2),'r');
%----------画初始点和最优点--------------------------------------
plot(x_initial(1),x_initial(2),'or')
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
function H_test1=H_test1(x)
x1=x(1);
x2=x(2);
h11=2/(x1^2+x2^2+2)^2-(8*x1^2)/(x1^2+x2^2+2)^3;
h12=-(8*x1*x2)/(x1^2+x2^2+2)^3;
h21=h12;
%连续函数的Hesse阵是对称阵
h22=2/(x1^2+x2^2+2)^2-(8*x2^2)/(x1^2+x2^2+2)^3;
H_test1=[h11,h12;h21,h22];
end