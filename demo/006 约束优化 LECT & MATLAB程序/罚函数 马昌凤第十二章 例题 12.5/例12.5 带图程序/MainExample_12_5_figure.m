function MainExample_12_5_figure()
close all;
clear all;
clc;
figure
hold on

kk=1:1:50;

x0=[1,3.5]';
[x,mu,lambda,output,fx_rec]=multphr_figure('f1','h1','g1','df1','dh1','dg1',x0);
n=length(fx_rec);
plot(kk(1:n),fx_rec(1:n),'color','b','Marker','.');
plot(kk(n),fx_rec(n),'color','b','Marker','p')

x0=[3,3]';
[x,mu,lambda,output,fx_rec]=multphr_figure('f1','h1','g1','df1','dh1','dg1',x0);
n=length(fx_rec);
plot(kk(1:n),fx_rec(1:n),'color','m','Marker','.');
plot(kk(n),fx_rec(n),'color','m','Marker','p')


x0=[1,1]';
[x,mu,lambda,output,fx_rec]=multphr_figure('f1','h1','g1','df1','dh1','dg1',x0);
n=length(fx_rec);
plot(kk(1:n),fx_rec(1:n),'color','r','Marker','.');
plot(kk(n),fx_rec(n),'color','r','Marker','p')


x0=[2,4]';
[x,mu,lambda,output,fx_rec]=multphr_figure('f1','h1','g1','df1','dh1','dg1',x0);
n=length(fx_rec);
plot(kk(1:n),fx_rec(1:n),'color','c','Marker','.');
plot(kk(n),fx_rec(n),'color','c','Marker','p')
end
%目标函数文件f1.m
function f=f1(x)
f=(x(1)-2.0)^2+(x(2)-1.0)^2;
end
%目标函数的梯度文件df1.m
function g=df1(x)
g=[2.0*(x(1)-2.0), 2.0*(x(2)-1.0)]';
end
%不等式约束函数文件g1.m
function gi=g1(x)
gi=-0.25*x(1)^2-x(2)^2+1;
end
%不等式约束（向量）函数的Jacobi矩阵（转置）文件dg1.m
function dgi=dg1(x)
dgi=[-0.5*x(1), -2.0*x(2)]';
end
%等式约束函数文件h1.m
function he=h1(x)
he=x(1)-2.0*x(2)+1.0;
end
%等式约束（向量）函数的Jacobi矩阵（转置）文件dh1.m
function dhe=dh1(x)
dhe=[1.0, -2.0]';
end
%%%%%%%%%%%%%%%%%% 增广拉格朗日函数%%%%%%%%%%%%%%%%%%%%%
function psi=mpsi(x,fun,hf,gf,dfun,dhf,dgf,mu,lambda,sigma)
f=feval(fun,x); he=feval(hf,x); gi=feval(gf,x);
l=length(he); m=length(gi);
psi=f; s1=0.0;
for i=1:l
psi=psi-he(i)*mu(i);
s1=s1+he(i)^2;
end
psi=psi+0.5*sigma*s1;
s2=0.0;
for i=1:m
s3=max(0.0, lambda(i) - sigma*gi(i));
s2=s2+s3^2-lambda(i)^2;
end
psi=psi+s2/(2.0*sigma);
end
%%%%%%%%%%%%%%%%%% 增广拉格朗日函数的梯度%%%%%%%%%%%%%%%%%%%%
function dpsi=dmpsi(x,fun,hf,gf,dfun,dhf,dgf,mu,lambda,sigma)
dpsi=feval(dfun,x);
he=feval(hf,x); gi=feval(gf,x);
dhe=feval(dhf,x); dgi=feval(dgf,x);
l=length(he); m=length(gi);
for i=1:l
  dpsi=dpsi+(sigma*he(i)-mu(i))*dhe(:,i);
end
for i=1:m
  dpsi=dpsi+(sigma*gi(i)-lambda(i))*dgi(:,i);
end
end
function [x,mu,lambda,output,fx_rec]=multphr_figure(fun,hf,gf,dfun,dhf,dgf,x0)
%功能: 用乘子法解一般约束问题: min f(x), s.t. h(x)=0, g(x)?=0
%输入: x0是初始点, fun, dfun分别是目标函数及其梯度；
% hf, dhf分别是等式约束（向量）函数及其Jacobi矩阵的转置；
% gf, dgf分别是不等式约束（向量）函数及其Jacobi矩阵的转置；
%输出: x是近似最优点，mu, lambda分别是相应于等式约束和不
% 等式约束的乘子向量; output是结构变量, 输出近似极小值f, 迭
% 代次数, 内迭代次数等

%======================
  %记录过程点的函数值
k_rec=1;
fx_rec=[];
fx_rec(1)=feval(fun,x0);
%======================

maxk=500; %最大迭代次数
sigma=2.0; %罚因子
eta=2.0; theta=0.8; %PHR算法中的实参数
k=0; ink=0; %k, ink分别是外迭代和内迭代次数
epsilon=1e-5; %终止误差值
x=x0; he=feval(hf,x); gi=feval(gf,x);
n=length(x); l=length(he); m=length(gi);
%选取乘子向量的初始值
mu=0.1*ones(l,1); lambda=0.1*ones(m,1);
btak=10; btaold=10; %用来检验终止条件的两个值
while(btak>epsilon && k<maxk)
%调用BFGS算法程序求解无约束子问题
  [x,ival,ik]=bfgs_phr('mpsi','dmpsi',x0,fun,hf,gf,dfun,dhf,dgf,mu,lambda,sigma);
  ink=ink+ik;
  he=feval(hf,x); gi=feval(gf,x);
  btak=0.0;
  for i=1:l
    btak=btak+he(i)^2; 
  end
  for i=1:m
    temp=min(gi(i),lambda(i)/sigma);
    btak=btak+temp^2;
  end
  btak=sqrt(btak);
  if btak>epsilon
    if(k>=2 && btak>theta*btaold)
       sigma=eta*sigma;
    end
%更新乘子向量
    for i=1:l, mu(i)=mu(i)-sigma*he(i); end
    for i=1:m
       lambda(i)=max(0.0,lambda(i)-sigma*gi(i));
    end
  end
  k=k+1;
  btaold=btak;
  x0=x;
  %======================
  %记录过程点的函数值
  k_rec=k_rec+1;
  fx_rec(k_rec)=feval(fun,x);
  %===================
end
f=feval(fun,x);
output.fval=f;
output.iter=k;
output.inner_iter=ink;
output.bta=btak;
end