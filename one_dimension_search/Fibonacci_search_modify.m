clc
clear;
%%  这是测试区
alpha_lower = 0;
alpha_upper = 1;
tolerance = 10^-4;
x_current = -0.5;
d_current = 1;
[alpha_star] = Fibonacci_search(@ftest1,alpha_lower,alpha_upper,tolerance,x_current,d_current);
%%  斐波那契（Fibonacci_search）搜索函数
function [alpha_star] = Fibonacci_search(f_test,alpha_lower,alpha_upper,tolerance,x_current,d_current)
    Fabonacci_series_upper = (alpha_upper - alpha_lower)/tolerance;
    Fabonacci_series = [1,2];
    n = 2;
    while(Fabonacci_series(n)<=Fabonacci_series_upper)
        n = n+1;
        Fabonacci_series(n) = Fabonacci_series(n-1)+Fabonacci_series(n-2);
    end
    %   k = 0时的左右点以及
    k = 0;
    reduction_rate_k = Fabonacci_series(n-1)/Fabonacci_series(n);
    alpha_upper_k = alpha_upper;
    alpha_lower_k = alpha_lower;
    L_k = (alpha_upper_k - alpha_lower_k);
    alpha_left_k = alpha_upper_k - L_k*reduction_rate_k;
    alpha_right_k = alpha_lower_k + L_k*reduction_rate_k;
    x_current_left_k = x_current + alpha_left_k*d_current;
    x_current_right_k = x_current + alpha_right_k*d_current;
    f_x_left_k = f_test(x_current_left_k);
    f_x_right_k = f_test(x_current_right_k);
    %   k>=1时开始迭代
    while(abs(alpha_right_k - alpha_left_k)>tolerance&&n>=4)
        n = n - 1;
        reduction_rate_k = Fabonacci_series(n-1)/Fabonacci_series(n);
        if(f_x_left_k >f_x_right_k)
            alpha_lower_k = alpha_left_k;
            alpha_left_k = alpha_right_k;
            L_k = abs(alpha_upper_k - alpha_lower_k);
            alpha_right_k = alpha_lower_k + L_k*reduction_rate_k;
        elseif(f_x_left_k <f_x_right_k)
            alpha_upper_k = alpha_right_k;
            alpha_right_k = alpha_left_k;
            L_k = abs(alpha_upper_k - alpha_lower_k);
            alpha_left_k = alpha_upper_k - L_k*reduction_rate_k;
        else
            alpha_lower_k = alpha_left_k;
            alpha_upper_k = alpha_right_k;
            L_k = abs(alpha_upper_k - alpha_lower_k);
            alpha_right_k = alpha_lower_k + L_k*reduction_rate_k;
            alpha_left_k = alpha_upper_k - L_k*reduction_rate_k;
        end
        x_current_left_k = x_current + alpha_left_k*d_current;
        x_current_right_k = x_current + alpha_right_k*d_current;
        f_x_left_k = f_test(x_current_left_k);
        f_x_right_k = f_test(x_current_right_k);
        k = k+1;
    end
    if(tolerance >10^8)
        disturbance_quantity = 10^-9;
    else
        disturbance_quantity = tolerance*0.1;
    end
    alpha_middle = (alpha_upper_k+alpha_lower_k)/2;
    alpha_right_k = alpha_middle + disturbance_quantity;
    alpha_left_k = alpha_middle - disturbance_quantity;
    x_current_left_k = x_current + alpha_left_k*d_current;
    x_current_right_k = x_current + alpha_right_k*d_current;
    f_x_right_k = f_test(x_current_right_k);
    f_x_left_k = f_test(x_current_left_k);
    if(f_x_left_k <f_x_right_k)
        alpha_upper_k = alpha_right_k;
    elseif(f_x_left_k >f_x_right_k)
        alpha_lower_k = alpha_left_k;
    else
        alpha_upper_k = alpha_rght_k;
        alpha_lower_k = alpha_left_k;
    end
    alpha_star = (alpha_upper_k + alpha_lower_k)/2;
end
%%  写一个测试函数
function ftest1 = ftest1(x)
    ftest1 = 2*x^2-x-1;
end



















