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
% ������ʽ
%[x_optimal,f_optimal,k,lamda_optimal,mu_optimal]=QP_acset_eqineq(H,p,q,Aeq,beq,Aineq,bineq,x_lower,x_upper)
%--------------------------------------------------------------------------
% �������˵��
%--------------------------------------------------------------------------
% H��͹���ι滮�����Hesse����
% p��͹���ι滮����һ����ϵ������
% q��͹���ι滮����ĳ�����
% Aeq�����Ե�ʽԼ����ϵ��������������ɳڱ������˹�������ϵ��
% beq�����Ե�ʽԼ����Ӧ���Ҷ�����
% Aineq:���Բ���ʽԼ����ϵ��������������ɳڱ������˹�������ϵ��,��>=��
% bineq:���Բ���ʽԼ����Ӧ���Ҷ�����,��>=����ʽ
% x_lower�������½�
% x_upper�������Ͻ�
%--------------------------------------------------------------------------
%�������
%--------------------------------------------------------------------------
%x_optimal�����ŵ�
%f_optimal��f_test��Ӧx_optimal�ĺ���ֵ
%k:������Ž���Ҫ�����Ĵ���
%lamda_optimal����Ӧ��ʽԼ�����������ճ�������ֵ
%mu_optimal����Ӧ����ʽԼ�����������ճ�������ֵ
%==========================================================================
global x_rec;
[m1,n1]=size(Aeq);
[m2,n2]=size(Aineq);
n=max(n1,n2);
%------------���ʼ���н�---------------------------------------------------
[x_initial,f_initial]=LP_simplex_phase1_general(Aeq,beq,Aineq,bineq,x_lower,x_upper);
%--------------------------------------------------------------------------

%--���������ϡ��½�ת�ɲ���ʽԼ������ԭ����ʽԼ����ɡ��µġ�һ�鲻��ʽԼ��-----
[Aineq,bineq]=LowerUpperToIneq(Aineq,bineq,x_lower,x_upper,n);
[m2,n2]=size(Aineq);
m=m1+m2;

%--------------------------------------------------------------------------
% �ҳ�x_initial���Ļ�����
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
    %  �⡰�ȼ����⡱�����d_active
    %----------------------------------------------------------------------
    g_k=H*x_k+p;
    A_active=A(acset,:);
    b_active=b(acset);
    if m_active>m1
        %--------�ҳ�����Լ���е������޹��ǲ��ֻ���Լ��-----------------------
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
    % �ж�����������Ӧ��������Լ�����Լ��ǻ���Լ����
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
%--------�ҳ�����Լ���е������޹��ǲ��ֻ���Լ��-----------------------
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
%   ���������ϡ��½�ת�ɲ���ʽԼ��
%=========================================================================
% n������ά��
%-----------------------------------------------------------------------
[m2,n2]=size(Aineq);
%--------------------------------------------------------------------------
%--------------���ϡ��½��Լ���벻��ʽԼ����ϳɡ��¡�����ʽԼ��-------------
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
% �����Թ滮�����η�-��һ�׶� ��� ���е�
%==========================================================================
function [x_initial,f_initial]=LP_simplex_phase1_general(Aeq,beq,Aineq,bineq,x_lower,x_upper)
%==========================================================================
%�������ø�ʽ��
%[x_initial,f_initial]=LP_simplex_phase1_general(Aeq,beq,Aineq,bineq,x_lower,x_upper)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%������ı��ʽԼ��Ϊ��>=��
% ���ϡ��½�תΪ����ʽԼ��
% �����ɳڱ��� n_with_slack=n+m2;
% [m1,n1]=size(Aeq); [m2.n2]=size(Aineq);n=max(n1,n2);m=m1+m2
% Aeq=[Aeq,zreos(m1,m2)];Aineq=[Aineq,-eye(m2,m2)];
% ת��ʹ���Ҷ�������Ԫ�طǸ�
% A=[Aeq;Aineq];A=[A,eye(m,m)];
% b=[beq;bineq]
% �����˹�����n_with_artificial=n_with_slack+m;_
% �����ε�һ�׶ε�Ŀ�꺯��ϵ����c_phase1=[n+m2��0��m1+m2��1];
% phase1 ���Ż����⣬ֻ�е�ʽԼ����û�в���ʽԼ��
%-----------------------------------------------------------------------
%�������˵��
%--------------------------------------------------------------------------
%Aeq�����Ե�ʽԼ����ϵ��������������ɳڱ������˹�������ϵ��
%beq�����Ե�ʽԼ����Ӧ���Ҷ�����
%Aineq:���Բ���ʽԼ����ϵ��������������ɳڱ������˹�������ϵ��
%bineq:���Բ���ʽԼ����Ӧ���Ҷ�����
%x_lower�������½�
%x_upper�������Ͻ�
%--------------------------------------------------------------------------
%�������
%--------------------------------------------------------------------------
%x_initial�����ŵ�
%f_initial��f_test��Ӧx_initial�ĺ���ֵ��Ϊ��
%==========================================================================
%==========================================================================
%������˵��
%==========================================================================
% ���ܣ���MATLAB�к���linprog�����м��е�һ�����н�
% ���м���
%  Aeq*x = beq
%  Aineq*x >= bineq
%  x_lower <= x <= x_upper
%--------------------------------------------------------------------------
%Ҫ������ı��ʽԼ��Ϊ��>=��
%--------------------------------------------------------------------------

%-----------------------------------------------------------------------
[m1,n1]=size(Aeq);
[m2,n2]=size(Aineq);
n=max(n1,n2);
%--------------���ϡ��½��Լ���벻��ʽԼ����ϳɡ��¡�����ʽԼ��-------------
[Aineq,bineq]=LowerUpperToIneq(Aineq,bineq,x_lower,x_upper,n);
%--------------------------------------------------------------------------
[m2,n2]=size(Aineq);
m=m1+m2;


%--------------------����ɳڱ�����Ҫ���ɳڱ�������Ǹ�Ҫ��-------------------
if (m1>0&&m2>0)
    Aeq=[Aeq,zeros(m1,m2)];
end
if (m2>0)
    Aineq=[Aineq,-eye(m2,m2)];
end
n_with_slack=n+m2;
%-------------------���Ҷ�������Ϊ�Ǹ���ʹ�˹���������Ǹ�Ҫ��---------------
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
%------------��ת����ĵ�ʽԼ������˹��������˹�����Ҫ��Ǹ�-----------------
m=m1+m2;
A=[Aeq;Aineq];
A=[A,eye(m,m)];
b=[beq;bineq];
n_with_artificial=n_with_slack+m;

%-----------��һ�׶����⿪ʼ------------------------------------------------
k=0;
base_set=(n_with_slack+1):n_with_artificial;
nonbase_set=setdiff(1:n_with_artificial,base_set);
B_inverse=eye(m,m);
N=A(:,nonbase_set);
c_phase1=[zeros(n_with_slack,1);ones(m,1)];
c_B=c_phase1(base_set);
c_N=c_phase1(nonbase_set);
lamda=(B_inverse')*c_B;

while(tolerance>0) %�����һ�׶�����ĵ����ε��� 
    %-----------���Ӧ�ڷǻ������ļ���������---------------------------------
    mu=c_N-(N')*lamda;
    %------------�жϵ�һ�׶������Ŀ�꺯���Ƿ��Ѿ��ﵽ����-------------------
    if(mu(:)>-tolerance)
        x_B=B_inverse*b;
        f_phase1=(c_B')*x_B;
        f_initial=(c_B')*x_B;
        break
    end
    %----------ȷ��������������ţ�ʹ��Ŀ�꺯�����ٵ����---------------------
    min_mu=0;
    for j=1:n_with_slack
        if(mu(j)<min_mu)
            min_mu=mu(j);
            index_entering_relative=j;
        end
    end
    index_entering=nonbase_set(index_entering_relative);
    %------��һ�׶����ⲻ������޽�������ֻ��ֱ��ȷ�����������Ӧ�����------
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
    %-------------------------����-----------------------------------------
    base_set(index_leaving_relative)=index_entering;
    nonbase_set(index_entering_relative)=index_leaving;
    %----------------------------------------------------------------------
    %-----------------------���µĻ�����������----------------------------
    k=k+1;
    B_inverse=B_inverse_update(B_inverse,vector_entering,index_leaving_relative,m);
    %-----------------------Ϊ��һ�ֵ�������--------------------------------
    N=A(:,nonbase_set);
    c_B=c_phase1(base_set);
    c_N=c_phase1(nonbase_set);
    lamda=(B_inverse')*c_B;
end
%-------------------�ж�ԭ�����Ƿ���ڻ������л�-----------------------------
if(f_phase1>tolerance)
    disp('The problem has no basic feasible solution');
    x_initial=NaN;
    f_initial=NaN;
    return
end

%----------���˹�������Ӧ����Ŵ�base_set�н�������--------------------------
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
%       �����ʼ���е�
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
%�������˵��
%--------------------------------------------------------------------------
%B_inversae��������B�������
%vector_entering����������
%r:����������B�е�������
%m��B_inversae��ά����������B_inversae�Ƿ���
%--------------------------------------------------------------------------
%�������˵��
%--------------------------------------------------------------------------
%B_inversae����������µĻ�����B�������
%--------------------------------------------------------------------------
vector_entering_B=B_inverse(r,:);
%B_inverse�е���һ�в��ܱ��滻
%--------------------------------------------------------------------------
%  r=1,�������m=1�����
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
%-----------------r=m,����ֻ����m>1�����-----------------------------------
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
%-------------r>1&r<m������ֻ����m>1�����----------------------------------
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
%    �������˵��
%====================================================================
% H��͹���ι滮�����HESSE����
% p��͹���ι滮�����һ����ϵ������
% q��Ŀ�꺯���ĳ�����
% Aeq�����Ե�ʽԼ����ϵ������
% beq�����Ե�ʽ���Ҷ�����
%=====================================================================
%      �������
% x_optimal,���Ž�����
% f_optimal,����ֵ
% k,������Ž�ĵ�������
% lamda_optimal����Ӧ�ڵ�ʽԼ�����������ճ��ӵ����Ž�
%=====================================================================
[m,n]=size(H);
%-------------���HESSE�������С����ֵ���ж��Ƿ���Ҫ����----------
epsilong=1e-10;
min_eigenvalue=min(eig(H));
if min_eigenvalue<-epsilong %�ж�����ֵ�Ƿ�С��0
    disp('the problem is not convex');
    x_optimal=NaN;
    f_optimal=NaN;
    k=NaN;
    lamda_optimal=NaN;
    return;
end
if abs(min_eigenvalue)<epsilong %�ж�����ֵ�Ƿ����0
    H=H+epsilong*eye(n,n);
end
%=================================================================
%    �����Լ����Aeq=[]����Hx=-p����x_optimal
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
%   ����QR�ֽ⽫ֻ����ʽԼ����͹���ι滮����ת������Լ����͹���ι滮����
%----------------------------------------------------------------------
[m,n]=size(Aeq);
k=1;
[Q,R]=qr(Aeq');%ע�⣺�ԡ�Aeq��ת�á���QR�ֽ�
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
%    ����LDLT�ֽⷨ��������������ճ���
%-----------------------------------------------------------------
[L,D]=ldl(Aeq*Aeq');
 lamda_optimal=(L')\(D\(L\(Aeq*(H*x_optimal+p))));
 %----------------------------------------------------------
end