function example_10_4_XinggaoLiu
close all;
clear all;
clc;

H=[1,-4,2,1;-4,16,-8,-4;2,-8,4,2;1,-4,2,1];
p=[-1;0;7;4];
q=5;
Aeq=[1,1,1,1];
beq=[4];
Aineq=[-1,-2,0,0];
bineq=[-3.5];
x_lower=zeros(4,1);
x_upper=[inf;0.5;inf;inf];

[x_optimal,f_optimal,k,lamda_optimal,mu_optimal]=QP_acset_eqineq(H,p,q,Aeq,beq,Aineq,bineq,x_lower,x_upper)

end



function [x_optimal,f_optimal,k,lamda_optimal,mu_optimal]=QP_acset_eqineq(H,p,q,Aeq,beq,Aineq,bineq,x_lower,x_upper)
%==========================================================================
% 调用形式
%[x_optimal,f_optimal,k,lamda_optimal,mu_optimal]=QP_acset_eqineq(H,p,q,Aeq,beq,Aineq,bineq,x_lower,x_upper)
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
% x_lower：变量下界
% x_upper：变量上界
%--------------------------------------------------------------------------
%输出参数
%--------------------------------------------------------------------------
%x_optimal：最优点
%f_optimal：f_test对应x_optimal的函数值
%k:获得最优解需要迭代的次数
%lamda_optimal：对应等式约束的拉格朗日乘子最优值
%mu_optimal：对应不等式约束的拉格朗日乘子最优值
%==========================================================================
global x_rec;
[m1,n1]=size(Aeq);
[m2,n2]=size(Aineq);
n=max(n1,n2);
%------------求初始可行解---------------------------------------------------
[x_initial,f_initial]=LP_simplex_phase1_general(Aeq,beq,Aineq,bineq,x_lower,x_upper);
%--------------------------------------------------------------------------

%--将变量的上、下界转成不等式约束，与原不等式约束组成“新的”一组不等式约束-----
[Aineq,bineq]=LowerUpperToIneq(Aineq,bineq,x_lower,x_upper,n);
[m2,n2]=size(Aineq);
m=m1+m2;

%--------------------------------------------------------------------------
% 找出x_initial处的积极集
%--------------------------------------------------------------------------
epsilong=1e-10;
A=[Aeq;Aineq];
b=[beq;bineq];
acset=1:m1;
m_active=m1;
k=0;
x_k=x_initial;

vector_temp=A((m1+1):m,:)*x_k-b((m1+1):m);
for i=1:m2
    if abs(vector_temp(i))<epsilong
        acset=[acset,m1+i];
        m_active=m_active+1;
    end
end
nonacset=setdiff(1:m,acset);
lamda_optimal=zeros(m,1);
while epsilong<1
    %----------------------------------------------------------------------
    %  解“等价问题”，获得d_active
    %----------------------------------------------------------------------
    g_k=H*x_k+p;
    A_active=A(acset,:);
    b_active=b(acset);
    if m_active>m1
        %--------找出积极约束中的线性无关那部分积极约束-----------------------
        [indexes_independent,number_independent]=Constraint_independent_eq(A_active,b_active,m_active,n);
        %------------------------------------------------------------------
        if number_independent<m_active
            m_active=number_independent;
            A_active=A_active(indexes_independent,:);
            acset=acset(indexes_independent);
            nonacset=setdiff(1:m,acset);
        end
    end
    [d_active,f_active,k_active,lamda_active]=QP_QR_eq(H,g_k,0,A_active,zeros(m_active,1));
    %----------------------------------------------------------------------
    % 判断收敛，并相应调整积极约束集以及非积极约束集
    %----------------------------------------------------------------------
    if abs(d_active)<epsilong
        if m_active>m1
            lamda_active_ineq=lamda_active((m1+1):m_active);
        else
            lamda_active_ineq=0;
        end
        if lamda_active_ineq(:)>-epsilong
            x_optimal=x_k;
            f_optimal=0.5*(x_optimal)'*(H*x_optimal+2*p)+q;
            lamda_optimal(acset)=lamda_active;
            mu_optimal=lamda_optimal((m1+1):m);
            lamda_optimal=lamda_optimal(1:m1);
            return;
        else
            min_lamda_active_ineq=lamda_active(m1+1);
            index_acset_relative=m1+1;
            for i=(m1+1):m_active
                lamda_temp=lamda_active(i);
                if lamda_temp<min_lamda_active_ineq
                    min_lamda_active_ineq=lamda_temp;
                    index_acset_relative=i;
                end
            end
            index_acset=acset(index_acset_relative);
            nonacset=union(nonacset,index_acset);
            acset=setdiff(acset,index_acset);
            m_active=m_active-1;
            k=k+1;
        end
    else
        min_alpha_k=1;
        m_nonactive=m-m_active;
        for i=1:m_nonactive
            index_nonacset=nonacset(i);
            scalar_temp=A(index_nonacset,:)*d_active;
            if scalar_temp<-epsilong
                alpha_k=(b(index_nonacset)-A(index_nonacset,:)*x_k)/scalar_temp;
                if alpha_k<min_alpha_k
                    min_alpha_k=alpha_k;
                    index_nonacset_relative=i;
                end
            end
        end
        alpha_k=min_alpha_k;
        x_k=x_k+alpha_k*d_active;
        k=k+1;

        if alpha_k<1
            index_nonacset=nonacset(index_nonacset_relative);
            acset=union(acset,index_nonacset);
            nonacset=setdiff(nonacset,index_nonacset);
            m_active=m_active+1;
        end
    end
end
%==========================================================================
end

function [indexes_independent,number_independent]=Constraint_independent_eq(A_active,b_active,m_active,n)
%--------找出积极约束中的线性无关那部分积极约束-----------------------
%[indexes_independent,number_independent]=Constraint_independent_eq(A_active,b_active,m_active,n);
%------------------------------------------------------------------
if m_active==0
    indexes_independent=[];
    number_independent=[];
else
    indexes_independent=[1];
    number_independent=1;
    if m_active<2
        return
    else
        Ab_matrix_current_0=[A_active(1,:),b_active(1)];
        for i=2:m_active
            Ab_row_add=[A_active(i,:),b_active(i)];
            Ab_matrix_current=[Ab_matrix_current_0;Ab_row_add];
            if rank(Ab_matrix_current)<number_independent+1
                Ab_matrix_current=Ab_matrix_current_0;
            else
                Ab_matrix_current_0=Ab_matrix_current;
                indexes_independent=[indexes_independent,i];
                number_independent=number_independent+1;
            end
        end
    end
    if number_independent>n
        disp('over-constrain');
    end
end   
end



function [Aineq,bineq]=LowerUpperToIneq(Aineq,bineq,x_lower,x_upper,n)
%=========================================================================
%   将变量的上、下界转成不等式约束
%=========================================================================
% n：变量维数
%-----------------------------------------------------------------------
[m2,n2]=size(Aineq);
%--------------------------------------------------------------------------
%--------------将上、下界的约束与不等式约束组合成“新”不等式约束-------------
unit_matrix=eye(n,n);
Aineq_add_x_lower=[];
bineq_add_x_lower=[];
for i=1:n
    if x_lower(i)~=-inf
        Aineq_add_x_lower=[Aineq_add_x_lower;unit_matrix(i,:)];
        bineq_add_x_lower=[bineq_add_x_lower;x_lower(i)];
        m2=m2+1;
    end
end
Aineq=[Aineq;Aineq_add_x_lower];
bineq=[bineq;bineq_add_x_lower];

Aineq_add_x_upper=[];
bineq_add_x_upper=[];
for i=1:n
    if x_upper(i)~=inf
        Aineq_add_x_upper=[Aineq_add_x_upper;-unit_matrix(i,:)];
        bineq_add_x_upper=[bineq_add_x_upper;-x_upper(i)];
        m2=m2+1;
    end
end
Aineq=[Aineq;Aineq_add_x_upper];
bineq=[bineq;bineq_add_x_upper];
end






%==========================================================================
% 用线性规划单纯形法-第一阶段 求解 可行点
%==========================================================================
function [x_initial,f_initial]=LP_simplex_phase1_general(Aeq,beq,Aineq,bineq,x_lower,x_upper)
%==========================================================================
%函数调用格式：
%[x_initial,f_initial]=LP_simplex_phase1_general(Aeq,beq,Aineq,bineq,x_lower,x_upper)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%将输入的表达式约束为“>=”
% 将上、下界转为不等式约束
% 增加松弛变量 n_with_slack=n+m2;
% [m1,n1]=size(Aeq); [m2.n2]=size(Aineq);n=max(n1,n2);m=m1+m2
% Aeq=[Aeq,zreos(m1,m2)];Aineq=[Aineq,-eye(m2,m2)];
% 转换使得右端向量的元素非负
% A=[Aeq;Aineq];A=[A,eye(m,m)];
% b=[beq;bineq]
% 增加人工变量n_with_artificial=n_with_slack+m;_
% 单纯形第一阶段的目标函数系数，c_phase1=[n+m2个0；m1+m2个1];
% phase1 的优化问题，只有等式约束，没有不等式约束
%-----------------------------------------------------------------------
%输入参数说明
%--------------------------------------------------------------------------
%Aeq：线性等式约束的系数矩阵，无需包含松弛变量和人工变量的系数
%beq：线性等式约束对应的右端向量
%Aineq:线性不等式约束的系数矩阵，无需包含松弛变量和人工变量的系数
%bineq:线性不等式约束对应的右端向量
%x_lower：变量下界
%x_upper：变量上界
%--------------------------------------------------------------------------
%输出参数
%--------------------------------------------------------------------------
%x_initial：最优点
%f_initial：f_test对应x_initial的函数值，为零
%==========================================================================
%==========================================================================
%主程序及说明
%==========================================================================
% 功能：用MATLAB中函数linprog求解可行集中的一个可行解
% 可行集：
%  Aeq*x = beq
%  Aineq*x >= bineq
%  x_lower <= x <= x_upper
%--------------------------------------------------------------------------
%要求输入的表达式约束为“>=”
%--------------------------------------------------------------------------

%-----------------------------------------------------------------------
[m1,n1]=size(Aeq);
[m2,n2]=size(Aineq);
n=max(n1,n2);
%--------------将上、下界的约束与不等式约束组合成“新”不等式约束-------------
[Aineq,bineq]=LowerUpperToIneq(Aineq,bineq,x_lower,x_upper,n);
%--------------------------------------------------------------------------
[m2,n2]=size(Aineq);
m=m1+m2;


%--------------------添加松弛变量，要求松弛变量满足非负要求-------------------
if (m1>0&&m2>0)
    Aeq=[Aeq,zeros(m1,m2)];
end
if (m2>0)
    Aineq=[Aineq,-eye(m2,m2)];
end
n_with_slack=n+m2;
%-------------------将右端向量化为非负，使人工变量满足非负要求---------------
tolerance=1e-14;
for i=1:m1
    if(beq(i)<-tolerance)
        beq(i)=-beq(i);
        Aeq(i,:)=-Aeq(i,:);
    end
end

for i=1:m2
    if(bineq(i)<-tolerance)
        bineq(i)=-bineq(i);
        Aineq(i,:)=-Aineq(i,:);
    end
end
%------------给转换后的等式约束添加人工变量，人工变量要求非负-----------------
m=m1+m2;
A=[Aeq;Aineq];
A=[A,eye(m,m)];
b=[beq;bineq];
n_with_artificial=n_with_slack+m;

%-----------第一阶段问题开始------------------------------------------------
k=0;
base_set=(n_with_slack+1):n_with_artificial;
nonbase_set=setdiff(1:n_with_artificial,base_set);
B_inverse=eye(m,m);
N=A(:,nonbase_set);
c_phase1=[zeros(n_with_slack,1);ones(m,1)];
c_B=c_phase1(base_set);
c_N=c_phase1(nonbase_set);
lamda=(B_inverse')*c_B;

while(tolerance>0) %进入第一阶段问题的单纯形迭代 
    %-----------求对应于非基变量的检验数向量---------------------------------
    mu=c_N-(N')*lamda;
    %------------判断第一阶段问题的目标函数是否已经达到最优-------------------
    if(mu(:)>-tolerance)
        x_B=B_inverse*b;
        f_phase1=(c_B')*x_B;
        f_initial=(c_B')*x_B;
        break
    end
    %----------确定进基变量的序号，使得目标函数减少的最多---------------------
    min_mu=0;
    for j=1:n_with_slack
        if(mu(j)<min_mu)
            min_mu=mu(j);
            index_entering_relative=j;
        end
    end
    index_entering=nonbase_set(index_entering_relative);
    %------第一阶段问题不会出现无界的情况，只需直接确定离基变量对应的序号------
    vector_entering=B_inverse*A(:,index_entering);
    x_B=B_inverse*b;
    delta=inf;
    for i=1:m
        if(vector_entering(i)>tolerance)
            delta_temp=x_B(i)/vector_entering(i);
            if(delta_temp<delta)
                delta=delta_temp;
                index_leaving_relative=i;
            end
        end
    end
    index_leaving=base_set(index_leaving_relative);
    %-------------------------换基-----------------------------------------
    base_set(index_leaving_relative)=index_entering;
    nonbase_set(index_entering_relative)=index_leaving;
    %----------------------------------------------------------------------
    %-----------------------求新的基矩阵的逆矩阵----------------------------
    k=k+1;
    B_inverse=B_inverse_update(B_inverse,vector_entering,index_leaving_relative,m);
    %-----------------------为下一轮迭代更新--------------------------------
    N=A(:,nonbase_set);
    c_B=c_phase1(base_set);
    c_N=c_phase1(nonbase_set);
    lamda=(B_inverse')*c_B;
end
%-------------------判断原问题是否存在基本可行基-----------------------------
if(f_phase1>tolerance)
    disp('The problem has no basic feasible solution');
    x_initial=NaN;
    f_initial=NaN;
    return
end

%----------将人工变量对应的序号从base_set中交换出来--------------------------
for i=1:m
    if(base_set(i)>n_with_slack)
        index_leaving_relative=i;
        index_leaving=base_set(index_leaving_relative);
        vector_swap=B_inverse(i,:)*N;
        for j=1:n_with_slack
            index_entering_relative=j;
            index_entering=nonbase_set(index_entering_relative);
            if((index_entering<=n_with_slack)&&(abs(vector_swap(j))>tolerance))
                break;
            end
        end
        base_set(index_leaving_relative)=index_ertering;
        nonbase_set(index_entering_relative)=index_leaving;
        %B_inverse=B_inverse_update(B_inverse,vector_entering,index_swap_leaving_relative,m);???
        B_inverse=B_inverse_update(B_inverse,vector_entering,index_leaving_relative,m);
    end
end
%----------------------------------------------------------------------------
%       计算初始可行点
%--------------------------------------------------------------------------
x_initial_0=zeros(n_with_slack,1);
x_B=B_inverse*b;
f_initial=(c_B')*x_B;
x_initial_0(base_set)=x_B;
x_initial=x_initial_0(1:n);
%---------------------------------------------------------------------
end


%==========================================================================
function B_inverse=B_inverse_update(B_inverse,vector_entering,r,m)
%--------------------------------------------------------------------------
%输入参数说明
%--------------------------------------------------------------------------
%B_inversae：基矩阵B的逆矩阵
%vector_entering：进基变量
%r:出基向量在B中的相对序号
%m：B_inversae的维数（阶数，B_inversae是方阵）
%--------------------------------------------------------------------------
%输出参数说明
%--------------------------------------------------------------------------
%B_inversae：换基后的新的基矩阵B的逆矩阵
%--------------------------------------------------------------------------
vector_entering_B=B_inverse(r,:);
%B_inverse中的这一行不能被替换
%--------------------------------------------------------------------------
%  r=1,这里包含m=1的情况
%--------------------------------------------------------------------------
if(r==1)
    %-----------------------i=1--------------------------------------------
    scalar_temp=vector_entering(r);
    for j=1:m
        B_inverse(r,j)=vector_entering_B(j)/scalar_temp;
    end
    %---------------------i+2~m--------------------------------------------
    if(m>1)
        for i=2:m
           scalar_temp=vector_entering(i)/vector_entering(r);
           for j=1:m
               B_inverse(i,j)=B_inverse(i,j)-scalar_temp*vector_entering_B(j);
           end
        end
    end
    return
end
%-----------------r=m,这里只包含m>1的情况-----------------------------------
if(r==m)
    %-------------i=1~m-1--------------------------------------------------
    for i=1:(m-1)
       scalar_temp=vector_entering(i)/vector_entering(r);
       for j=1:m
          B_inverse(i,j)=B_inverse(i,j)-scalar_temp*vector_entering_B(j);
       end
    end 
    %--------------------i=m-----------------------------------------------
    scalar_temp=vector_entering(r);
    for j=1:m
       B_inverse(r,j)=vector_entering_B(j)/scalar_temp;
    end
    return;
end
%-------------r>1&r<m，这里只包含m>1的情况----------------------------------
%--------------i=1~r-1-----------------------------------------------------
 for i=1:(r-1)
    scalar_temp=vector_entering(i)/vector_entering(r);
    for j=1:m
       B_inverse(i,j)=B_inverse(i,j)-scalar_temp*vector_entering_B(j);
    end
 end
 %-------------------i=r---------------------------------------------------
 scalar_temp=vector_entering(r);
 for j=1:m
    B_inverse(r,j)=vector_entering_B(j)/scalar_temp;
 end
 %---------------------i=r+1~m---------------------------------------------
 for i=(r+1):m
     scalar_temp=vector_entering(i)/vector_entering(r);
    for j=1:m
       B_inverse(i,j)=B_inverse(i,j)-scalar_temp*vector_entering_B(j);
    end
 end
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