function [alpha_star,x_next,f_next,k]=Dichotomous_search(f_test,x_current,d_current,alpha_lower,alpha_upper,tolerance)
%==========================================================================
%     �������ø�ʽ
%[alpha_star,x_next,f_next,k]=Dichotomous_search(@f_test,x_current,d_current,alpha_lower,alpha_upper,tolerance)
%--------------------------------------------------------------------------
%   �������˵��
%--------------------------------------------------------------------------
%f_test��Ŀ�꺯��
%x_current��x�������ռ��еĵ�ǰ�㣨��ȷ����
%d_current��f_test��x_current���½�����������ȷ����
%alpha_lower����x_current����������d_current�õ�һ��������������󣬸õ�������ĳ�ʼ�½�
%alpha_upper:��x_current����������d_current�õ�һ��������������󣬸õ�������ĳ�ʼ�Ͻ�
%tolerance�����������Ҫ�󳤶ȣ����Ȳ���С���Ŷ���
%--------------------------------------------------------------------------
%�������
%--------------------------------------------------------------------------
%alpha_star����ɶԷ�������õ��Ĳ���
%x_next��=,x_current+alpha_star*d_current
%f_next��f_test�ڵ�x_next����ֵ
%k����ɶԷ���������ĵ�������
%==========================================================================
%==========================================================================
% ������˵��
%--------------------------------------------------------------------------
%�Է�������ÿ�������ж���Ҫ���������²�������Ӧ������x�㣩���亯��ֵ
%k���Է��������еĴ���
%alpha_left_k������k�ζԷ�����ʱ�����е㸺�Ŷ���Ĳ�����f_alpha_left_kΪ��Ӧ�ĺ���ֵ
%alpha_right_k������k�ζԷ�����ʱ�����е㸺�Ŷ���Ĳ�����f_alpha_right_kΪ��Ӧ�ĺ���ֵ
%alpha_middle_k������[alpha_lower_k��alpha_upper_k]���е㣬f_alpha_middle_kΪ��Ӧ�ĺ���ֵ
%alpha_lower_k:����k�ζԷ��������õĲ����������˵㣬f_alpha_lower_kΪ��Ӧ�ĺ���ֵ
%alpha_upper_k:����k�ζԷ��������õĲ����������˵㣬f_alpha_upper_kΪ��Ӧ�ĺ���ֵ
%disturbance_quantity:��alph_middle��΢С�Ŷ�����ԽС��ÿ����ȥ������Խ�ӽ��ڵ�ǰ�����1/2
%  disturbance_quantity <= (1/2)tolerance�����򣬶Է�����������������
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%΢С�Ŷ����ļ��㣬�Ŷ���Ĭ��Ϊ1e-9����tolerance�ǳ�Сʱ���Ŷ�������Ϊtolerance��1/2����
%--------------------------------------------------------------------------
if(tolerance>=1e-8)
    disturbance_quantity=1e-9;
else
    disturbance_quantity=0.1*tolerance;
end
%--------------------------------------------------------------------------
% k=0ʱ
%--------------------------------------------------------------------------
k=0;
alpha_lower_k=alpha_lower;
alpha_upper_k=alpha_upper;
%--------------------------------------------------------------------------
% k>=1ʱ
%--------------------------------------------------------------------------
while (abs(alpha_upper_k-alpha_lower_k)>tolerance)
    alpha_middle_k=0.5*(alpha_upper_k+alpha_lower_k);
    alpha_left_k=alpha_middle_k-disturbance_quantity;
    alpha_right_k=alpha_middle_k+disturbance_quantity;
    x_alpha_left_k=x_current+alpha_left_k*d_current;
    x_alpha_right_k=x_current+alpha_right_k*d_current;
    f_alpha_left_k=f_test(x_alpha_left_k);
    f_alpha_right_k=f_test(x_alpha_right_k);
    if(f_alpha_left_k<f_alpha_right_k)
        alpha_upper_k=alpha_right_k;
    elseif(f_alpha_left_k>f_alpha_right_k)
        alpha_lower_k=alpha_left_k;
    else
        alpha_lower_k=alpha_left_k;
        alpha_upper_k=alpha_right_k;
    end
    k=k+1;
end
%--------------------------------------------------------------------------
% ȡ����[alpha_lower_k��alpha_upper_k]���е���Ϊ��Ѳ���
%--------------------------------------------------------------------------
alpha_star=0.5*(alpha_lower_k+alpha_upper_k);
x_next=x_current+alpha_star*d_current;
f_next=f_test(x_next);
end



