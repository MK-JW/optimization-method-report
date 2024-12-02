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
% ������ʽ
%[x_optimal,f_optimal,k,lamda_optimal,mu_optimal]=QP_pdi_infeasible_eqineq(H,p,q,Aeq,beq,Aineq,bineq,z_z_initial,tolerance)
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
% z_initial��ԭ����ż��ĳ�ʼ����ֵ
% tolerance����ż����ľ���Ҫ��
%--------------------------------------------------------------------------
%�������
%--------------------------------------------------------------------------
%x_optimal�����ŵ�
%f_optimal����Ӧx_optimal�ĺ���ֵ
%k:������Ž���Ҫ�����Ĵ���
%lamda_optimal����Ӧ��ʽԼ�����������ճ�������ֵ
%mu_optimal����Ӧ����ʽԼ�����������ճ�������ֵ
%==========================================================================
global x_rec;

[m1,n1]=size(Aeq);
[m2,n2]=size(Aineq);
n=max(n1,n2);
%------------��ʼ�����㷨����-----------------------------------------
k=0;
reduce_factor_tau=m2/(m2+10*sqrt(m2));
reduce_factor_alpha=1-1e-6;
x_initial=z_initial(1:n);
y_initial=z_initial(n+1:n+m2);
lamda_initial=z_initial(n+m2+1:n+m1+m2);
mu_initial=z_initial(n+m1+m2+1:n+m1+2*m2);
%--------------------------------------------------------------------------
%�ж�y_initial��mu_initial��Ԫ���Ƿ�Ϊ��;��<=0,�����Ԫ��=1
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

%================��¼������======================
   x_rec(1,:)=x_k;
%================================================

duality_gap_k=max(1,y_k'*mu_k);
tau_k=duality_gap_k/m2;
while duality_gap_k>tolerance
    %--------���³ͷ�����-----------
    tau_next=tau_k*reduce_factor_tau;
    %--------�����µķ���----------------------------------------------------
    vector_mu=y_k-Aineq*x_k+bineq;
    %----------------------------------------------------------------------
    %-------------���С��ޣ�������ʽԼ���ֱ���-----------------------------
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
   %-------------���С��޵�ʽԼ���ֱ���------------------------------------ 
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
%  ��ԭ����Ͷ�ż��������µĲ���
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
%     ������һ����
%--------------------------------------------------------------------------
   x_k=x_k+alpha_k*delta_x;
   y_k=y_k+alpha_k*delta_y;
   %-------------���С��޵�ʽԼ���ֱ���------------------------------------
   if isempty(lamda_k)
       lamda_k=[];%���޵�ʽԼ��ʱ
   else
       lamda_k=lamda_k+alpha_k*delta_lamda;
   end
   %-----------------------------------------------------------------------
   mu_k=mu_k+alpha_k*delta_mu;
%--------------------------------------------------------------------------
%   �����µĶ�ż���ֵ
%--------------------------------------------------------------------------
   duality_gap_k=y_k'*mu_k;
   tau_k=duality_gap_k/m2;
   k=k+1;
   %================��¼������======================
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