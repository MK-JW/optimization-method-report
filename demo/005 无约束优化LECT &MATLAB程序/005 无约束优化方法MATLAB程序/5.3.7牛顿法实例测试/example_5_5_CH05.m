function example_5_5_CH05
close;
clear;
clc;
%-----------------------------------------------
%���˸�-Ӧ�����Ż�������MATLABʵ��CH05-5.2
%��5.5  Ŀ�꺯��f(x)=x1^2+x2^2+x1*x2+2;
%  x0=(1;-4);  tol=1e-6
%------------------------------------------------
x_initial=[1;-4];
tolerance=1e-6;
[x_optimal,f_optimal,k]=Newton(@f_test2,@g_test2,@H_test2,x_initial,tolerance)
t=-5:0.1:5;
[X1,X2]=meshgrid(t);
Y=X1.^2+X2.^2+X1.*X2+2;
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
function H_test2=H_test2(x)
x1=x(1);
x2=x(2);
h11=2;
h12=1;
h21=h12;
%����������Hesse���ǶԳ���
h22=2;
H_test2=[h11,h12;h21,h22];
end