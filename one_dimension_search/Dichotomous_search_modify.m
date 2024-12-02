clc
clear
%%  测试一下计算的步长
alpha_lower = 0;
alpha_upper = 1;
tolerance = 10^-4;
x_current = -0.5;
d_current = 1;
[alpha_star] = Dichotomous_search(@f_test1,alpha_lower,alpha_upper,tolerance,x_current,d_current);
%%  对分搜索法（Dichotomous_search）计算步长
function [alpha_star]=Dichotomous_search(f_test,alpha_lower,alpha_upper,tolerance,x_current,d_current)
    if(tolerance>=10^-8)
        disturbance_quantity = 10^-9;
    else
        disturbance_quantity = tolerance*0.1;
    end
    %   k = 0 时的情况
    k = 0;
    alpha_lower_k = alpha_lower;
    alpha_upper_k = alpha_upper;
    alpha_left_k = (alpha_lower_k+alpha_upper_k)/2 - disturbance_quantity;
    alpha_right_k = (alpha_lower_k+alpha_upper_k)/2 + disturbance_quantity;
    x_current_left_k = x_current+alpha_left_k*d_current;
    x_current_right_k = x_current+alpha_right_k*d_current;
    f_alpha_left_k = f_test(x_current_left_k);
    f_alpha_right_k = f_test(x_current_right_k);
    %   k >= 1 时的情况
    while(abs(alpha_upper_k - alpha_lower_k)>tolerance)
        if(f_alpha_left_k<f_alpha_right_k)
            alpha_upper_k = alpha_right_k;
        elseif(f_alpha_left_k>f_alpha_right_k)
            alpha_lower_k = alpha_left_k;
        else
            alpha_upper_k = alpha_right_k;
            alpha_lower_k = alpha_left_k;
        end
        alpha_left_k = (alpha_lower_k+alpha_upper_k)/2 - disturbance_quantity;
        alpha_right_k = (alpha_lower_k+alpha_upper_k)/2 + disturbance_quantity;
        x_current_left_k = x_current+alpha_left_k*d_current;
        x_current_right_k = x_current+alpha_right_k*d_current;
        f_alpha_left_k = f_test(x_current_left_k);
        f_alpha_right_k = f_test(x_current_right_k);
        k = k+1;
    end
    alpha_star = (alpha_upper_k+alpha_lower_k)/2;
end
%%  定义一个测试函数f_test用来测试计算
function f_test1 = f_test1(x)
    f_test1 = 2*x^2 - x -1;
end
