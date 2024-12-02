function [alpha_star,x_next,f_next,k]=Trisection_search(f_test,x_current,d_current,alpha_lower,alpha_upper,tolerance)
%==========================================================================
%     �������ø�ʽ
%[alpha_star,x_next,f_next,k]=Trisection_search(@f_test,x_current,d_current,alpha_lower,alpha_upper,tolerance)
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
%alpha_star��������ȷ�������õ��Ĳ���
%x_next��=,x_current+alpha_star*d_current
%f_next��f_test�ڵ�x_next����ֵ
%k��������ȷ���������ĵ�������
%==========================================================================
%==========================================================================
% ������˵��
%--------------------------------------------------------------------------
%���ȷ�������ÿ�������ж���Ҫ���������²�������Ӧ������x�㣩���亯��ֵ
%k�����ȷ��������еĴ���
%alpha_left_k������k�����ȷ�����ʱ������[alpha_lower_k��alpha_upper_k]���ȷֵ�����ߵĵ㣬f_alpha_left_kΪ��Ӧ�ĺ���ֵ
%alpha_right_k������k�����ȷ�����ʱ������[alpha_lower_k��alpha_upper_k]���ȷֵ����ұߵĵ㣬f_alpha_right_kΪ��Ӧ�ĺ���ֵ
%alpha_middle_k������[alpha_lower_k��alpha_upper_k]���ȷֵ��м�㣬f_alpha_middle_kΪ��Ӧ�ĺ���ֵ
%alpha_lower_k:����k�ζԷ��������õĲ����������˵㣬f_alpha_lower_kΪ��Ӧ�ĺ���ֵ
%alpha_upper_k:����k�ζԷ��������õĲ����������˵㣬f_alpha_upper_kΪ��Ӧ�ĺ���ֵ
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%k=0ʱ����Ҫ���κ���ֵ����
%--------------------------------------------------------------------------
k=0;
alpha_lower_k=alpha_lower;
alpha_upper_k=alpha_upper;
alpha_left_k=alpha_lower_k+(1/4)*(alpha_upper_k-alpha_lower_k);
x_alpha_left_k=x_current+alpha_left_k*d_current;
f_alpha_left_k=f_test(x_alpha_left_k);
alpha_middle_k=alpha_lower_k+(1/2)*(alpha_upper_k-alpha_lower_k);
x_alpha_middle_k=x_current+alpha_middle_k*d_current;
f_alpha_middle_k=f_test(x_alpha_middle_k);
alpha_right_k=alpha_lower_k+(3/4)*(alpha_upper_k-alpha_lower_k);
x_alpha_right_k=x_current+alpha_right_k*d_current;
f_alpha_right_k=f_test(x_alpha_right_k);
%--------------------------------------------------------------------------
% k>=1ʱ
%--------------------------------------------------------------------------
while (abs(alpha_upper_k-alpha_lower_k)>tolerance)
    k=k+1;
    if (f_alpha_left_k<=f_alpha_middle_k)&&(f_alpha_left_k<=f_alpha_right_k)
    alpha_upper_k=alpha_middle_k;
    alpha_middle_k=alpha_left_k;
    f_alpha_middle_k=f_alpha_left_k;
    alpha_left_k=alpha_lower_k+(1/4)*(alpha_upper_k-alpha_lower_k);
    x_alpha_left_k=x_current+alpha_left_k*d_current;
    f_alpha_left_k=f_test(x_alpha_left_k);
    alpha_right_k=alpha_lower_k+(3/4)*(alpha_upper_k-alpha_lower_k);
    x_alpha_right_k=x_current+alpha_right_k*d_current;
    f_alpha_right_k=f_test(x_alpha_right_k);
    elseif(f_alpha_middle_k<=f_alpha_left_k)&&(f_alpha_middle_k<=f_alpha_right_k)
        alpha_lower_k=alpha_left_k;
        alpha_upper_k=alpha_right_k;
        alpha_left_k=alpha_lower_k+(1/4)*(alpha_upper_k-alpha_lower_k);
        x_alpha_left_k=x_current+alpha_left_k*d_current;
        f_alpha_left_k=f_test(x_alpha_left_k);
        alpha_right_k=alpha_lower_k+(3/4)*(alpha_upper_k-alpha_lower_k);
        x_alpha_right_k=x_current+alpha_right_k*d_current;
        f_alpha_right_k=f_test(x_alpha_right_k);
    else
        alpha_lower_k=alpha_middle_k;
        alpha_middle_k=alpha_right_k;
        f_alpha_middle_k=f_alpha_right_k;
        alpha_left_k=alpha_lower_k+(1/4)*(alpha_upper_k-alpha_lower_k);
        x_alpha_left_k=x_current+alpha_left_k*d_current;
        f_alpha_left_k=f_test(x_alpha_left_k);
        alpha_right_k=alpha_lower_k+(3/4)*(alpha_upper_k-alpha_lower_k);
        x_alpha_right_k=x_current+alpha_right_k*d_current;
        f_alpha_right_k=f_test(x_alpha_right_k);
    end
end
%--------------------------------------------------------------------------
% ȡ����[alpha_lower_k��alpha_upper_k]���е���Ϊ��Ѳ���
%--------------------------------------------------------------------------
alpha_star=0.5*(alpha_lower_k+alpha_upper_k);
x_next=x_current+alpha_star*d_current;
f_next=f_test(x_next);
end



