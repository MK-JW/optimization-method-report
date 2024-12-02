function example_11_3_XinggaoLiu
close all;
clear all;
clc;
global x_rec;
H=[4,-1,0,0;-1,2,0,0;0,0,2,0;0,0,0,8];
p=[-4;0;-2;0];
q=11;
Aeq=[1,2,-1,1];
beq=[0];
Aineq=[-1,1,-2,1;2,-2,3,3];
bineq=[-7;9];

x_initial=[2;-1;3;3];
[m1,n1]=size(Aeq);
[m2,n2]=size(Aineq);
z_initial=[x_initial;ones(m1+2*m2,1)];
tolerance=1e-4;
[x_optimal,f_optimal,k,lamda_optimal,mu_optimal]=QP_pdi_infeasible_eqineq(H,p,q,Aeq,beq,Aineq,bineq,z_initial,tolerance)
x_record=x_rec
end



function [x_optimal,f_optimal,k,lamda_optimal,mu_optimal]=QP_pdi_infeasible_eqineq(H,p,q,Aeq,beq,Aineq,bineq,z_initial,tolerance)
%==========================================================================
% 调用形式
%[x_optimal,f_optimal,k,lamda_optimal,mu_optimal]=QP_pdi_infeasible_eqineq(H,p,q,Aeq,beq,Aineq,bineq,z_z_initial,tolerance)
%--------------------------------------------------------------------------
% 输入参数说明
%--------------------------------------------------------------------------
% H：凸二次规划问题的Hesse矩阵
% p：凸二次规划问题一次项系数向量
% q：凸二次规划问题的常数项
% Aeq：线性等式约束的系数矩阵，无需包含松弛变量和人工变量的系数
% beq：线性等式约束对应的右端向量
% Aineq:线性不等式约束的系数矩阵，无需包含松弛变量和人工变量的系数,“>=”
% bineq:线性不等式约束对应的右端向量,“>=”形式
% z_initial：原—对偶解的初始估计值
% tolerance：对偶间隔的精度要求
%--------------------------------------------------------------------------
%输出参数
%--------------------------------------------------------------------------
%x_optimal：最优点
%f_optimal：对应x_optimal的函数值
%k:获得最优解需要迭代的次数
%lamda_optimal：对应等式约束的拉格朗日乘子最优值
%mu_optimal：对应不等式约束的拉格朗日乘子最优值
%==========================================================================
global x_rec;

[m1,n1]=size(Aeq);
[m2,n2]=size(Aineq);
n=max(n1,n2);
%------------开始进入算法处理-----------------------------------------
k=0;
reduce_factor_tau=m2/(m2+10*sqrt(m2));
reduce_factor_alpha=1-1e-6;
x_initial=z_initial(1:n);
y_initial=z_initial(n+1:n+m2);
lamda_initial=z_initial(n+m2+1:n+m1+m2);
mu_initial=z_initial(n+m1+m2+1:n+m1+2*m2);
%--------------------------------------------------------------------------
%判断y_initial和mu_initial的元素是否为正;若<=0,令非正元素=1
%--------------------------------------------------------------------------
for i=1:m2
    if y_initial(i)<=0
        y_initial(i)=1;
    end
    if mu_initial(i)<=0
        mu_initial(i)=1;
    end
end
x_k=x_initial;
y_k=y_initial;
lamda_k=lamda_initial;
mu_k=mu_initial;

%================记录迭代点======================
   x_rec(1,:)=x_k;
%================================================

duality_gap_k=max(1,y_k'*mu_k);
tau_k=duality_gap_k/m2;
while duality_gap_k>tolerance
    %--------更新惩罚因子-----------
    tau_next=tau_k*reduce_factor_tau;
    %--------计算新的方向----------------------------------------------------
    vector_mu=y_k-Aineq*x_k+bineq;
    %----------------------------------------------------------------------
    %-------------分有、无（不）等式约束分别处理-----------------------------
    Aeq_lamda_k=Aeq'*lamda_k;
    Aineq_mu_k=Aineq'*mu_k;
    if isempty(Aeq_lamda_k)
        Aeq_lamda_k=zeros(n,1);
    end
    if isempty(Aineq_mu_k)
        Aineq_mu_k=zeros(n,1);
    end
    %----------------------------------------------------------------------
    vector_x=H*x_k+p-Aeq_lamda_k-Aineq_mu_k;
    %------------------------------------------
    y_reciprocal=1./y_k;
    vector_y=tau_next*y_reciprocal-mu_k;
    vector_temp1=mu_k.*y_reciprocal;
    vector_temp2=vector_temp1.*vector_mu;
    Y_inverse_Mu=diag(vector_temp1);
    N=H+Aineq'*Y_inverse_Mu*Aineq;
    p_k=vector_x-Aineq'*(vector_y+vector_temp2);
   %-------------分有、无等式约束分别处理------------------------------------ 
   if isempty(Aeq)
       vector_lamda=0;
   else
       vector_lamda=Aeq*x_k-beq;
   end
   %----------------------------------------------------------------------- 
    [delta_x,f_delta_x,k,delta_lamda]=QP_QR_eq(N,p_k,0,Aeq,-vector_lamda);
    delta_y=Aineq*delta_x-vector_mu;
    delta_mu=vector_y-vector_temp1.*delta_y;
%--------------------------------------------------------------------------
%  对原问题和对偶问题计算新的步长
%--------------------------------------------------------------------------
   alpha_primal=1;
   alpha_dual=1;
   for i=1:m2
       if delta_y(i)<0
           alpha_temp=-y_k(i)/delta_y(i);
           if alpha_temp<alpha_primal
               alpha_primal=alpha_temp;
           end
       end
       if delta_mu(i)<0
           alpha_temp=-mu_k(i)/delta_mu(i);
           if alpha_temp<alpha_dual
               alpha_dual=alpha_temp;
           end
       end
   end
   alpha_primal=reduce_factor_alpha*alpha_primal;
   alpha_dual=reduce_factor_alpha*alpha_dual;
   alpha_k=min(alpha_primal,alpha_dual);
%--------------------------------------------------------------------------
%     计算下一个点
%--------------------------------------------------------------------------
   x_k=x_k+alpha_k*delta_x;
   y_k=y_k+alpha_k*delta_y;
   %-------------分有、无等式约束分别处理------------------------------------
   if isempty(lamda_k)
       lamda_k=[];%当无等式约束时
   else
       lamda_k=lamda_k+alpha_k*delta_lamda;
   end
   %-----------------------------------------------------------------------
   mu_k=mu_k+alpha_k*delta_mu;
%--------------------------------------------------------------------------
%   计算新的对偶间隔值
%--------------------------------------------------------------------------
   duality_gap_k=y_k'*mu_k;
   tau_k=duality_gap_k/m2;
   k=k+1;
   %================记录迭代点======================
   x_rec(k,:)=x_k;
   %================================================
end
x_optimal=x_k;
y_optimal=y_k;
lamda_optimal=lamda_k;
mu_optimal=mu_k;
f_optimal=0.5*x_optimal'*(H*x_optimal+2*p)+q;
%==========================================================================
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