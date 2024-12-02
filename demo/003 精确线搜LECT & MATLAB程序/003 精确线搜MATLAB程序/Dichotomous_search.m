function [alpha_star,x_next,f_next,k]=Dichotomous_search(f_test,x_current,d_current,alpha_lower,alpha_upper,tolerance)
%==========================================================================
%     函数调用格式
%[alpha_star,x_next,f_next,k]=Dichotomous_search(@f_test,x_current,d_current,alpha_lower,alpha_upper,tolerance)
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
%alpha_star：完成对分搜索后得到的步长
%x_next：=,x_current+alpha_star*d_current
%f_next：f_test在点x_next处的值
%k：完成对分搜索所需的迭代次数
%==========================================================================
%==========================================================================
% 主程序及说明
%--------------------------------------------------------------------------
%对分搜索在每轮搜索中都需要计算两个新步长（对应的两个x点）及其函数值
%k：对分搜索进行的次数
%alpha_left_k：经过k次对分搜索时，对中点负扰动后的步长，f_alpha_left_k为相应的函数值
%alpha_right_k：经过k次对分搜索时，对中点负扰动后的步长，f_alpha_right_k为相应的函数值
%alpha_middle_k：区间[alpha_lower_k，alpha_upper_k]的中点，f_alpha_middle_k为相应的函数值
%alpha_lower_k:经过k次对分搜索后获得的步长区间的左端点，f_alpha_lower_k为对应的函数值
%alpha_upper_k:经过k次对分搜索后获得的步长区间的左端点，f_alpha_upper_k为对应的函数值
%disturbance_quantity:对alph_middle的微小扰动量，越小则每次消去的区间越接近于当前区间的1/2
%  disturbance_quantity <= (1/2)tolerance，否则，对分搜索方法不收敛。
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%微小扰动量的计算，扰动量默认为1e-9，当tolerance非常小时将扰动量设置为tolerance的1/2以下
%--------------------------------------------------------------------------
if(tolerance>=1e-8)
    disturbance_quantity=1e-9;
else
    disturbance_quantity=0.1*tolerance;
end
%--------------------------------------------------------------------------
% k=0时
%--------------------------------------------------------------------------
k=0;
alpha_lower_k=alpha_lower;
alpha_upper_k=alpha_upper;
%--------------------------------------------------------------------------
% k>=1时
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
% 取区间[alpha_lower_k，alpha_upper_k]的中点作为最佳步长
%--------------------------------------------------------------------------
alpha_star=0.5*(alpha_lower_k+alpha_upper_k);
x_next=x_current+alpha_star*d_current;
f_next=f_test(x_next);
end



