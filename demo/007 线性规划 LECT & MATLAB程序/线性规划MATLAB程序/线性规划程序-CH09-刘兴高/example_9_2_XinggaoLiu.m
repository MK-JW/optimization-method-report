function example_9_2_XinggaoLiu
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
x_initial=[0.612;0.9268;2.0349];
z_initial=[x_initial;1.8015;1.9651;0.6120;0.9269;2.0349;1;5.75;2;1;1];
min_max=0;
tolerance=1e-4;
[x_optimal,y_optimal,f_optimal,k,lamda_optimal,mu_optimal]=LP_pdi_feasible_general(c,Aeq,beq,Aineq,bineq,z_initial,min_max,tolerance)
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

function [x_optimal,y_optimal,f_optimal,k,lamda_optimal,mu_optimal]=LP_pdi_feasible_general(c,Aeq,beq,Aineq,bineq,z_initial,min_max,tolerance)
%=============�������======================
% c:     Ŀ�꺯��ϵ������
% Aeq��  ��ʽԼ��ϵ������
% beq��  ��ʽԼ�����Ҷ�����
% Aineq������ʽԼ��ϵ������
% bineq������ʽԼ�����Ҷ�����
% z_initial: ԭ-��ż��ĳ�ʼֵ
% min_max: min��ȡ1��max��ȡ0��Ĭ��Ϊmin��
% toleranc:��ż����ľ���Ҫ��
%===========================================
%================�������˵��====================
% x_optimal�����߱��������Ž�
% y_optimal����Ӧ�ڲ���ʽԼ�����ɳڱ��������Ž�
% f_optimal��Ŀ�꺯������ֵ
% k��������Ž�ĵ�������
% lamda_optimal����Ӧ�ڵ�ʽԼ�����������ճ�������ֵ
% mu_optimal����Ӧ�ڲ���ʽԼ�����������ճ�������ֵ
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
x_k=x_initial;
y_k=y_initial;
lamda_k=lamda_initial;
mu_k=mu_initial;
duality_gap=max(1,(y_k')*mu_k);
tau_k=duality_gap/m2;
%-----��¼������-------
j=k+1;
x_rec(:,j)=x_k;
%-----------------------
while duality_gap>tolerance
    %====���·�����====
    tau_next=tau_k*reduce_factor_tau;
    %=====�����µķ���=========
    y_reciprocal=(1./y_k);
    vector_y=tau_next*y_reciprocal-mu_k;
    vector_temp1=mu_k.*y_reciprocal;
    Y_inverse_Mu=diag(vector_temp1);
    H=(Aineq')*Y_inverse_Mu*Aineq;
    p_k=-(Aineq')*vector_y;
    [delta_x,f_QP,exitflag,output,lambda]=quadprog(H,p_k,[],[],Aeq,zeros(m1,1),[],[]);
    delta_lamda=-lambda.eqlin;
    delta_y=Aineq*delta_x;
    delta_mu=vector_y-vector_temp1.*delta_y;
    %==��ԭ����Ͷ�ż����ֱ�����²���========
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
    %====������һ����=============================
    x_k=x_k+alpha_primal*delta_x;
    y_k=y_k+alpha_primal*delta_y;
    lamda_k=lamda_k+alpha_dual*delta_lamda;
    mu_k=mu_k+alpha_dual*delta_mu;
    %=====�����µĶ�ż���========================
    duality_gap=(y_k')*mu_k;
    tau_k=duality_gap/m2;
    k=k+1;
    %-----��¼������-------
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