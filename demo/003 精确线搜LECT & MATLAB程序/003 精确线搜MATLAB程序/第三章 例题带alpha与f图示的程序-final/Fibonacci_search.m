function [alpha_star,x_next,f_next,k]=Fibonacci_search(f_test,x_current,d_current,alpha_lower,alpha_upper,tolerance)
%==========================================================================
%     �������ø�ʽ
%[alpha_star,x_next,f_next,k]=Fibonacci_search(@f_test,x_current,d_current,alpha_lower,alpha_upper,tolerance)
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
%alpha_star�����Fibonacci������õ��Ĳ���
%x_next��=,x_current+alpha_star*d_current
%f_next��f_test�ڵ�x_next����ֵ
%k�����Fibonacci��������ĵ�������
%==========================================================================
%==========================================================================
% ������˵��
%--------------------------------------------------------------------------
%Fibonacci������ÿ�������ж���Ҫ���������²�������Ӧ������x�㣩���亯��ֵ
%k��Fibonacci�������еĴ���
%alpha_left_k������k��Fibonacci����ʱ������[alpha_lower_k��alpha_upper_k]Fibonacci����˵ĵ㣬f_alpha_left_kΪ��Ӧ�ĺ���ֵ
%alpha_right_k������k��Fibonacci����ʱ������[alpha_lower_k��alpha_upper_k]Fibonacci���Ҷ˵ĵ㣬f_alpha_right_kΪ��Ӧ�ĺ���ֵ
%length_k������k��Fibonacci���������Ĳ�������ĳ���
%n:ʹ�õ���Fibonacci���������ֵ���±�
%n-1:��Ҫ�������ܴ��������һ��Ϊ�Է�����
%k=n-1�������Ҫ���Է���������ΪFibonacci����F1/F2=1/2��
%���յõ��Ĳ������䳤�ȿ��ܻ��Դ���tolerance���������ټ������жԷ�������
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Fibonacci_series����ֵ�ļ��㼰nֵ����
%--------------------------------------------------------------------------
%Fibonacciʹ��Fibonacci_series��F1��ʼ��ֻ�����Fibonacci_series(1)~~(n)
%Fibonacci_series_upper=(alpha_upper-alpha_lower)/tolerance
%------------------------------------------------------------
Fibonacci_series_upper=(alpha_upper-alpha_lower)/tolerance;
Fibonacci_series=[1,2];
n=2;
while(Fibonacci_series(n)<=Fibonacci_series_upper)
    n=n+1;
    Fibonacci_series(n)=Fibonacci_series(n-1)+Fibonacci_series(n-2);
end
%--------------------------------------------------------------------------

%==============plot phi(alpha) ==================
N=100;N1=N+1;
d_a=(alpha_upper-alpha_lower)/N;
a=alpha_lower:d_a:alpha_upper;
f_a=ones(N1,1);
for i=1:N1
    x_a=x_current+a(i)*d_current;
    f_a(i,1)=f_test(x_a);
end
plot(a,f_a(:,1),'k');
hold on
%===========================================  

%================��(0,f(x_current))��(1,f(x_current+d_current===============
plot(0,f_test(x_current),'ob');
hold on
plot(1,f_test(x_current+d_current),'og');
xlabel('alpha');ylabel('f(x+alpha*d)');
%========================================================================
%--------------------------------------------------------------------------
%k=0ʱ����Ҫ��2�κ���ֵ����
%--------------------------------------------------------------------------
k=0;
alpha_lower_k=alpha_lower;
alpha_upper_k=alpha_upper;
length_k=(Fibonacci_series(n-1)/Fibonacci_series(n))*(alpha_upper_k-alpha_lower_k);
alpha_left_k=alpha_upper_k-length_k;
x_alpha_left_k=x_current+alpha_left_k*d_current;
f_alpha_left_k=f_test(x_alpha_left_k);
alpha_right_k=alpha_lower_k+length_k;
x_alpha_right_k=x_current+alpha_right_k*d_current;
f_alpha_right_k=f_test(x_alpha_right_k);

%=============��ͼ=========================================================    
   plot(alpha_left_k,f_alpha_left_k,'ob');
   plot(alpha_right_k,f_alpha_right_k,'og'); 
%=========================================================================

%--------------------------------------------------------------------------
% k=1~~n-2ʱ
%--------------------------------------------------------------------------
while (k<n-2)
    k=k+1;
    if (f_alpha_left_k<=f_alpha_right_k)
      alpha_upper_k=alpha_right_k;
      length_k=(Fibonacci_series(n-k-1)/Fibonacci_series(n-k))*(alpha_upper_k-alpha_lower_k);
      alpha_right_k=alpha_left_k;
      f_alpha_right_k=f_alpha_left_k;
      alpha_left_k=alpha_upper_k-length_k;
      x_alpha_left_k=x_current+alpha_left_k*d_current;
      f_alpha_left_k=f_test(x_alpha_left_k);
    else
        alpha_lower_k=alpha_left_k;
        length_k=(Fibonacci_series(n-k-1)/Fibonacci_series(n-k))*(alpha_upper_k-alpha_lower_k);
        alpha_left_k=alpha_right_k;
        f_alpha_left_k=f_alpha_right_k;
        alpha_right_k=alpha_lower_k+length_k;
        x_alpha_right_k=x_current+alpha_right_k*d_current;
        f_alpha_right_k=f_test(x_alpha_right_k);
    end
%=============��ͼ=========================================================    
   plot(alpha_left_k,f_alpha_left_k,'ob');
   plot(alpha_right_k,f_alpha_right_k,'og'); 
%=========================================================================
end
%--------------------------------------------------------------------------
% ���������whileѭ��k=n-2
%--------------------------------------------------------------------------
% k=n-1ʱ��alpha_left_k=alpha_right_k����̽�������1�㣬��ˣ����öԷַ���1������
%--------------------------------------------------------------------------
k=k+1;%k=n-1
alpha_middle_k=(1/2)*(alpha_lower_k+alpha_upper_k);
disturbance_quantity=0.1*tolerance;
alpha_left_k=alpha_middle_k-disturbance_quantity;
x_alpha_left_k=x_current+alpha_left_k*d_current;
f_alpha_left_k=f_test(x_alpha_left_k);
alpha_right_k=alpha_middle_k+disturbance_quantity;
x_alpha_right_k=x_current+alpha_right_k*d_current;
f_alpha_right_k=f_test(x_alpha_right_k);

%=============��ͼ=========================================================    
   plot(alpha_left_k,f_alpha_left_k,'ob');
   plot(alpha_right_k,f_alpha_right_k,'og'); 
%=========================================================================

if(f_alpha_left_k<f_alpha_right_k)
     alpha_upper_k=alpha_right_k;
elseif(f_alpha_left_k>f_alpha_right_k)
     alpha_lower_k=alpha_left_k;
else
     alpha_lower_k=alpha_left_k;
     alpha_upper_k=alpha_right_k;
end
%--------------------------------------------------------------------------
% ȡ����[alpha_lower_k��alpha_upper_k]���е���Ϊ��Ѳ���
%--------------------------------------------------------------------------
alpha_star=0.5*(alpha_lower_k+alpha_upper_k);
x_next=x_current+alpha_star*d_current;
f_next=f_test(x_next);
%=============��ͼ=========================================================    
   plot(alpha_star,f_next,'*r');
%=========================================================================
end



