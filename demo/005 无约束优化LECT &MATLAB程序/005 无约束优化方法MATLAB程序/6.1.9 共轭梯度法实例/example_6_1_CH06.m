function example_6_1_CH06
close;
clear;
clc;
%-----------------------------------------------
%���˸�-Ӧ�����Ż�������MATLABʵ��CH06-6.1
%��6.1  Ŀ�꺯��f(x)=-0.8*exp(-x1^2-4*x2^2);
%  x0=(1;1);  tol=1e-6
%------------------------------------------------
x_initial=[1;1];
tolerance=1e-6;
[x_optimal,f_optimal,k]=Conjugate_gradient_DY(@f_test1,@g_test1,x_initial,tolerance)
t=-3:0.1:3;
[X1,X2]=meshgrid(t);
Y=-0.8*exp(-X1.^2-4*X2.^2);
%---------------------------------------------------------------------
%��ȡ��¼��������ļ�����������
%---------------------------------------------------------------------
ex=importdata('testdata.txt');
data=ex.data;
[m,n]=size(data);
K=data(1:m,1)+1;%��������
X=data(1:m,2:n-1);%������
F=data(1:m,n);%������ĺ���ֵ
%-----------���ȸ��ߣ�ֵΪF�ģ�---------------
contour(X1,X2,Y,F)
%-------------------------------------------
%----------��������켣---------------------
hold on
plot(X(1:m,1),X(1:m,2),'r');
%----------����ʼ������ŵ�--------------------------------------
plot(x_initial(1),x_initial(2),'or')
plot(x_optimal(1),x_optimal(2),'*r')
end
function f_test1=f_test1(x)
x1=x(1);
x2=x(2);
f_test1=-0.8*exp(-x1^2-4*x2^2);
end
function g_test1=g_test1(x)
x1=x(1);
x2=x(2);
g1=(8*x1)/(5*exp(x1^2+4*x2^2));
g2=(32*x2)/(5*exp(x1^2+4*x2^2));
g_test1=[g1;g2];
end