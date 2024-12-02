function example_9_8_XinggaoLiu
close all;
clear all;
clc;
global x_rec
c=[1;1;5];
Aeq=[];
beq=[];
Aineq=[-3,-2,-1/4;0,0,-1];
bineq=[-6;-4];
Aineq=[Aineq;eye(3,3)];
bineq=[bineq;zeros(3,1)];
x_initial=[2.5;2.5;3];
[m1,n1]=size(Aeq);
[m2,n2]=size(Aineq);
z_initial=[x_initial;ones(m1+2*m2,1)];
min_max=0;
tolerance=1e-4;
[x_optimal,y_optimal,f_optimal,k,lamda_optimal,mu_optimal]=LP_pdi_pc_general(c,Aeq,beq,Aineq,bineq,z_initial,min_max,tolerance)
p1x=[0;5/3;2;0;0];p1y=[2.5;0;0;3;2.5];p1z=[4;4;0;0;4];
hold on
plot3(p1x,p1y,p1z,'k');
p2x=[0;5/3;0;0];p2y=[0;0;2.5;0];p2z=[4;4;4;4];
plot3(p2x,p2y,p2z,'k');
xlabel('x_1');ylabel('x_2');zlabel('x_3');
plot3(x_rec(1,:),x_rec(2,:),x_rec(3,:),'-o')
hold on
plot3(x_optimal(1),x_optimal(2),x_optimal(3),'*r')
end

function [x_optimal,y_optimal,f_optimal,k,lamda_optimal,mu_optimal]=LP_pdi_pc_general(c,Aeq,beq,Aineq,bineq,z_initial,min_max,tolerance)
%=============输入参数======================
% c:     目标函数系数向量
% Aeq：  等式约束系数矩阵
% beq：  等式约束的右端向量
% Aineq：不等式约束系数矩阵
% bineq：不等式约束的右端向量
% z_initial: 原-对偶解的初始值
% min_max: min型取1；max型取0；默认为min型
% toleranc:对偶间隔的精度要求
%===========================================
%================输出参数说明====================
% x_optimal：决策变量的最优解
% y_optimal：对应于不等式约束的松弛变量的最优解
% f_optimal：目标函数最优值
% k：获得最优解的迭代次数
% lamda_optimal：对应于等式约束的拉格朗日乘子最优值
% mu_optimal：对应于不等式约束的拉格朗日乘子最优值
%===================================================
global x_rec
[m1,n1]=size(Aeq);
[m2,n2]=size(Aineq);
n=max(n1,n2);
if min_max==0
    c=-c;
end

k=0;
reduce_factor_tau=m2/(m2+10*sqrt(m2));
reduce_factor_alpha=1-1e-6;

x_initial=z_initial(1:n);
y_initial=z_initial(n+1:n+m2);
lamda_initial=z_initial(n+m2+1:n+m1+m2);
mu_initial=z_initial(n+m1+m2+1:n+m1+m2+m2);
%======判断y_initial和mu_initial的元素是否为正，否则，令非正元素为1 ======
for i=1:m2
    if y_initial(i)<0
        y_initial(i)=1;
    end
    if mu_initial(i)<0
        mu_initial(i)=1;
    end
end
%=======================================================================
x_k=x_initial;
y_k=y_initial;
lamda_k=lamda_initial;
mu_k=mu_initial;
duality_gap=max(1,(y_k')*mu_k);
tau_k=duality_gap/m2;
%-----记录迭代点-------
j=k+1;
x_rec(:,j)=x_k;
%-----------------------
while duality_gap>tolerance   
    %=====计算新的方向=========
    vector_mu=y_k-Aineq*x_k+bineq;
    if isempty(Aeq)
        vector_x=c-(Aineq')*mu_k;
    else
        vector_x=c-(Aeq')*lamda_k-(Aineq')*mu_k;
    end
    y_reciprocal=(1./y_k);
    
    vector_temp1=mu_k.*y_reciprocal;
    vector_temp2=vector_temp1.*vector_mu;
    Y_inverse_Mu=diag(vector_temp1);
    H=(Aineq')*Y_inverse_Mu*Aineq;
    p_k_pre=-(Aineq')*(vector_temp2-mu_k)+vector_x;
    if isempty(Aeq)
        vector_lamda=[]
    else
        vector_lamda=Aeq*x_k-beq;
    end
    [delta_x_pre,f_QP,exitflag,output,lambda]=quadprog(H,p_k_pre,[],[],Aeq,-vector_lamda,[],[]);
    delta_lamda_pre=-lambda.eqlin;
    delta_y_pre=Aineq*delta_x_pre-vector_mu;
    delta_mu_pre=-mu_k-Y_inverse_Mu*delta_y_pre;
    %==对原问题和对偶问题分别计算新步长========
    alpha_primal_pre=1;
    alpha_dual_pre=1;
    for i=1:m2
        if delta_y_pre(i)<0
            alpha_temp=-(y_k(i)/delta_y_pre(i));
            if alpha_temp<alpha_primal_pre
                alpha_primal_pre=alpha_temp;
            end
        end
        if delta_mu_pre(i)<0
            alpha_temp=-(mu_k(i)/delta_mu_pre(i));
            if alpha_temp<alpha_dual_pre
                alpha_dual_pre=alpha_temp;
            end
        end
    end
    alpha_primal_pre=reduce_factor_alpha*alpha_primal_pre;
    alpha_dual_pre=reduce_factor_alpha*alpha_dual_pre;
    %====计算中心参数和缩减因子==============================
    tau_k=mu_k'*y_k/m2;
    tau_k_pre=[(mu_k+alpha_dual_pre*delta_mu_pre)'*(y_k+alpha_primal_pre*delta_y_pre)]/m2;
    sigma_k=(tau_k_pre/tau_k)^3;
    %====更新罚因子====
    tau_next=tau_k*sigma_k;
    %===计算校正方向=========================================
    vector_y=tau_next*y_reciprocal-(y_reciprocal.*delta_y_pre).*delta_mu_pre;
    H
    p_k_cor=-Aineq'*vector_y
    Aeq
   if isempty(Aeq)
       b_cor=[]
    else
        b_cor=zeros(m1,1)
    end
    [delta_x_cor,f_QP,exitflag,output,lambda]=quadprog(H,p_k_cor,[],[],Aeq,b_cor,[],[]);
    delta_lamda_cor=-lambda.eqlin;
    delta_y_cor=Aineq*delta_x_cor;
    delta_mu_cor=vector_y-Y_inverse_Mu*delta_y_cor;
    %=====计算搜索方向==================================
    delta_x=delta_x_pre+delta_x_cor;
    delta_y=delta_y_pre+delta_y_cor;
    delta_lamda=delta_lamda_pre+delta_lamda_cor;
    delta_mu=delta_mu_pre+delta_mu_cor;
    %=====计算步长======================================
    alpha_primal=1;
    alpha_dual=1;
    for i=1:m2
        if delta_y(i)<0
            alpha_temp=-(y_k(i)/delta_y(i));
            if alpha_temp<alpha_primal
                alpha_primal=alpha_temp;
            end
        end
        if delta_mu(i)<0
            alpha_temp=-(mu_k(i)/delta_mu(i));
            if alpha_temp<alpha_dual
                alpha_dual=alpha_temp;
            end
        end
    end
    alpha_primal=reduce_factor_alpha*alpha_primal;
    alpha_dual=reduce_factor_alpha*alpha_dual;    
    %====计算下一个点=============================
    x_k=x_k+alpha_primal*delta_x;
    y_k=y_k+alpha_primal*delta_y;
    lamda_k=lamda_k+alpha_dual*delta_lamda;
    mu_k=mu_k+alpha_dual*delta_mu;
    %=====计算新的对偶间隔========================
    duality_gap=(y_k')*mu_k;
    %tau_k=duality_gap/m2;
    k=k+1;
    %-----记录迭代点-------
    j=k+1;
    x_rec(:,j)=x_k;
    %----------------------
end
x_optimal=x_k;
y_optimal=y_k;
lamda_optimal=lamda_k;
mu_optimal=mu_k;
if min_max==0
    f_optimal=-(c')*x_k;
else
    f_optimal=(c')*x_k;
end
end