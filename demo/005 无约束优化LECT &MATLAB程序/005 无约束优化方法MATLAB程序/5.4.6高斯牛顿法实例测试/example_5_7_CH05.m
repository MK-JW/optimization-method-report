function example_5_7_CH05
close;
clear;
clc;
%-----------------------------------------------
%刘兴高-应用最优化方法及MATLAB实现CH05-5.7
%例5.7  目标函数F(x)=(x1^2+x2^2-1)^2+(x1+x2-2)^2;
%  [f1;f2]=[x1^2+x2^2-1;x1+x2-2]
%  x0=(2;2);  tol=1e-6
%------------------------------------------------
x_initial=[2;2];
%x_initial=[0.7487;0.7487];
tolerance=1e-6;
[x_optimal,F_optimal,k]=Guass_Newton(@F_test1,@G_test1,@J_test1,x_initial,tolerance)
t=-4:0.1:4;
[X1,X2]=meshgrid(t);
Y=(X1.^2+X2.^2-1).^2+(X1+X2-2).^2;

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
for i=1:m
plot(X(i,1),X(i,2),'*')
end
%----------画初始点和最优点--------------------------------------
plot(x_initial(1),x_initial(2),'or')
plot(x_optimal(1),x_optimal(2),'*G')
end
function F_test1=F_test1(x)
x1=x(1);
x2=x(2);
F_test1=(x1^2+x2^2-1)^2+(x1+x2-2)^2;
end
function G_test1=G_test1(x)
x1=x(1);
x2=x(2);
G1=4*x1*(x1^2+x2^2-1)+2*(x1+x2-2);
G2=4*x2*(x1^2+x2^2-1)+2*(x1+x2-2);
G_test1=[G1;G2];
end
function J_test1=J_test1(x)
x1=x(1);
x2=x(2);
J11=2*x1;
J12=2*x2;
J21=1;
J22=1;
J_test1=[J11,J12;J21,J22];
end