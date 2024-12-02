function [alpha_star,x_next,f_next,k]=Trisection_search(f_test,x_current,d_current,alpha_lower,alpha_upper,tolerance)
%==========================================================================
%     函数调用格式
%[alpha_star,x_next,f_next,k]=Trisection_search(@f_test,x_current,d_current,alpha_lower,alpha_upper,tolerance)
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
%alpha_star：完成三等分搜索后得到的步长
%x_next：=,x_current+alpha_star*d_current
%f_next：f_test在点x_next处的值
%k：完成三等分搜索所需的迭代次数
%==========================================================================
%==========================================================================
% 主程序及说明
%--------------------------------------------------------------------------
%三等分搜索在每轮搜索中都需要计算两个新步长（对应的两个x点）及其函数值
%k：三等分搜索进行的次数
%alpha_left_k：经过k次三等分搜索时，区间[alpha_lower_k，alpha_upper_k]三等分的最左边的点，f_alpha_left_k为相应的函数值
%alpha_right_k：经过k次三等分搜索时，区间[alpha_lower_k，alpha_upper_k]三等分的最右边的点，f_alpha_right_k为相应的函数值
%alpha_middle_k：区间[alpha_lower_k，alpha_upper_k]三等分的中间点，f_alpha_middle_k为相应的函数值
%alpha_lower_k:经过k次对分搜索后获得的步长区间的左端点，f_alpha_lower_k为对应的函数值
%alpha_upper_k:经过k次对分搜索后获得的步长区间的左端点，f_alpha_upper_k为对应的函数值
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
xlabel('alpha');ylabel('f(alpha)');
%========================================================================
%--------------------------------------------------------------------------
%k=0时，需要三次函数值估计
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
%=============画图=========================================================    
   plot(alpha_left_k,f_alpha_left_k,'ob');
    plot(alpha_middle_k,f_alpha_middle_k,'om')
   plot(alpha_right_k,f_alpha_right_k,'og'); 
%=========================================================================



% k>=1时
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
    %=============画图=========================================================    
   plot(alpha_left_k,f_alpha_left_k,'ob');
    plot(alpha_middle_k,f_alpha_middle_k,'om')
   plot(alpha_right_k,f_alpha_right_k,'og'); 
%=========================================================================
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



