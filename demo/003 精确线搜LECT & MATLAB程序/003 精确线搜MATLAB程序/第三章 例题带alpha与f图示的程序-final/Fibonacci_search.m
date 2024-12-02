function [alpha_star,x_next,f_next,k]=Fibonacci_search(f_test,x_current,d_current,alpha_lower,alpha_upper,tolerance)
%==========================================================================
%     函数调用格式
%[alpha_star,x_next,f_next,k]=Fibonacci_search(@f_test,x_current,d_current,alpha_lower,alpha_upper,tolerance)
%--------------------------------------------------------------------------
%   输入参数说明
%--------------------------------------------------------------------------
%f_test：目标函数
%x_current：x在向量空间中的当前点（已确定）
%d_current：f_test在x_current的下降搜索方向（已确定）
%alpha_lower：从x_current出发，沿着d_current得到一个步长单峰区间后，该单峰区间的初始下界
%alpha_upper:从x_current出发，沿着d_current得到一个步长单峰区间后，该单峰区间的初始上界
%tolerance：最终区间的要求长度，精度不能小于扰动量
%--------------------------------------------------------------------------
%输出参数
%--------------------------------------------------------------------------
%alpha_star：完成Fibonacci搜索后得到的步长
%x_next：=,x_current+alpha_star*d_current
%f_next：f_test在点x_next处的值
%k：完成Fibonacci搜索所需的迭代次数
%==========================================================================
%==========================================================================
% 主程序及说明
%--------------------------------------------------------------------------
%Fibonacci搜索在每轮搜索中都需要计算两个新步长（对应的两个x点）及其函数值
%k：Fibonacci搜索进行的次数
%alpha_left_k：经过k次Fibonacci搜索时，区间[alpha_lower_k，alpha_upper_k]Fibonacci的左端的点，f_alpha_left_k为相应的函数值
%alpha_right_k：经过k次Fibonacci搜索时，区间[alpha_lower_k，alpha_upper_k]Fibonacci的右端的点，f_alpha_right_k为相应的函数值
%length_k：经过k次Fibonacci搜索缩减的不定区间的长度
%n:使用到的Fibonacci序列中最大值的下标
%n-1:需要搜索的总次数，最后一次为对分搜索
%k=n-1的情况需要做对分搜索，因为Fibonacci序列F1/F2=1/2。
%最终得到的步长区间长度可能会略大于tolerance，本程序不再继续进行对分搜索。
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Fibonacci_series序列值的计算及n值估计
%--------------------------------------------------------------------------
%Fibonacci使用Fibonacci_series从F1开始，只需计算Fibonacci_series(1)~~(n)
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

%================画(0,f(x_current))和(1,f(x_current+d_current===============
plot(0,f_test(x_current),'ob');
hold on
plot(1,f_test(x_current+d_current),'og');
xlabel('alpha');ylabel('f(x+alpha*d)');
%========================================================================
%--------------------------------------------------------------------------
%k=0时，需要做2次函数值估计
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

%=============画图=========================================================    
   plot(alpha_left_k,f_alpha_left_k,'ob');
   plot(alpha_right_k,f_alpha_right_k,'og'); 
%=========================================================================

%--------------------------------------------------------------------------
% k=1~~n-2时
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
%=============画图=========================================================    
   plot(alpha_left_k,f_alpha_left_k,'ob');
   plot(alpha_right_k,f_alpha_right_k,'og'); 
%=========================================================================
end
%--------------------------------------------------------------------------
% 经过上面的while循环k=n-2
%--------------------------------------------------------------------------
% k=n-1时，alpha_left_k=alpha_right_k，试探步长变成1点，因此，采用对分法做1次搜索
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

%=============画图=========================================================    
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
% 取区间[alpha_lower_k，alpha_upper_k]的中点作为最佳步长
%--------------------------------------------------------------------------
alpha_star=0.5*(alpha_lower_k+alpha_upper_k);
x_next=x_current+alpha_star*d_current;
f_next=f_test(x_next);
%=============画图=========================================================    
   plot(alpha_star,f_next,'*r');
%=========================================================================
end



