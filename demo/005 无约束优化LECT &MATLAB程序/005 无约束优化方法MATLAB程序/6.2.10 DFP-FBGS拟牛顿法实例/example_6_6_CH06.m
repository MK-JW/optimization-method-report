function example_6_6_CH06
close;
clear;
clc;
%-----------------------------------------------
%刘兴高-应用最优化方法及MATLAB实现CH06-6.6
%例6.6  目标函数f(x)=x1^2+x2^2+x1*x2+2;
%  x0=(1;-4);  tol=1e-6
%------------------------------------------------
x_initial=[1;-4];
tolerance=1e-6;
[x_optimal,f_optimal,k]=BFGS_Wolfe(@f_test2,@g_test2,x_initial,tolerance)
t=-4.5:0.1:4.5;
[X1,X2]=meshgrid(t);
Y=X1.^2+X2.^2+X1.*X2+2;
%---------------------------------------------------------------------
%读取记录迭代点的文件，并画出点
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