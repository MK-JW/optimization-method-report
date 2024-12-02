function example_10_2_XinggaoLiu
close all;
clear all;
clc;
H=eye(3);
p=[2;1;-1];
q=3;
Aeq=[3,-1,2;-2,4,0;-4,3,8];
beq=[7;12;10];
[x_optimal,f_optimal,k,lamda_optimal]=QP_QR_eq(H,p,q,Aeq,beq)
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
    