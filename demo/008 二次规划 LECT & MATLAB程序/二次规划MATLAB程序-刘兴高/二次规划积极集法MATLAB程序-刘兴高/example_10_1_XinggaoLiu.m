function example_10_1_XinggaoLiu
close all;
clear all;
clc;
H=[4,1;1,2];
p=[-1;-1];
q=2;
Aeq=[1,1];
beq=[1];
[x_optimal,f_optimal,k,lamda_optimal]=QP_QR_eq(H,p,q,Aeq,beq)
%==========画图=====================
t1=-1.5:0.05:2;
t2=-1.5:0.05:2.5;
[X,Y]=meshgrid(t1,t2);
Z=2*X.*X+Y.*Y+X.*Y-X-Y+2;
contour(X,Y,Z)
hold on
v=[f_optimal,f_optimal];
contour(X,Y,Z,v)
plot(x_optimal(1),x_optimal(2),'*r')
xc=[-1.5;2];yc=[2.5;-1];
plot(xc,yc,'k')
colorbar
end
function [x_optimal,f_optimal,k,lamda_optimal]=QP_QR_eq(H,p,q,Aeq,beq)
%====================================================================
%    输入参数说明
%====================================================================
% H：凸二次规划问题的HESSE矩阵
% p：凸二次规划问题的一次项系数向量
% q：目标函数的常数项
% Aeq：线性等式约束的系数矩阵
% beq：线性等式的右端向量
%=====================================================================
%      输出参数
% x_optimal,最优解向量
% f_optimal,最优值
% k,获得最优解的迭代次数
% lamda_optimal，对应于等式约束的拉格朗日乘子的最优解
%=====================================================================
[m,n]=size(H);
%-------------获得HESSE矩阵的最小特征值，判断是否需要修正----------
epsilong=1e-10;
min_eigenvalue=min(eig(H));
if min_eigenvalue<-epsilong %判断特征值是否小于0
    disp('the problem is not convex');
    x_optimal=NaN;
    f_optimal=NaN;
    k=NaN;
    lamda_optimal=NaN;
    return;
end
if abs(min_eigenvalue)<epsilong %判断特征值是否等于0
    H=H+epsilong*eye(n,n);
end
%=================================================================
%    如果无约束即Aeq=[]，解Hx=-p，得x_optimal
%-----------------------------------------------------------------
if isempty(Aeq)
    lamda_optimal=[];
    [L,D]=ldl(H);
    x_optimal=(L')\(D\(L\(-p)));
    k=1;
    f_optimal=0.5*((x_optimal')*H*x_optimal)+p'*x_optimal+q;
    return
end
%=====================================================================
%   采用QR分解将只含等式约束的凸二次规划问题转化成无约束的凸二次规划问题
%----------------------------------------------------------------------
[m,n]=size(Aeq);
k=1;
[Q,R]=qr(Aeq');%注意：对“Aeq的转置”做QR分解
Q1=Q(:,1:m);
R=R(1:m,:);
x_part1=Q1*inv(R)'*beq;
if m==n
    x_optimal=x_part1;
    f_optimal=0.5*((x_optimal')*H*x_optimal)+p'*x_optimal+q;
else
    Q2=Q(:,m+1:n);
    H_wave=(Q2')*H*(Q2);
    p_wave=(Q2')*(H*x_part1+p);
    [L,D]=ldl(H_wave);
    phi=(L')\(D\(L\(-p_wave)));
    x_part2=Q2*phi;
    x_optimal=x_part1+x_part2;
    f_optimal=0.5*((x_optimal')*H*x_optimal)+p'*x_optimal+q;
end
%==================================================================
%    利用LDLT分解法求解最优拉格朗日乘子
%-----------------------------------------------------------------
[L,D]=ldl(Aeq*Aeq');
 lamda_optimal=(L')\(D\(L\(Aeq*(H*x_optimal+p))));
 %----------------------------------------------------------
end
    