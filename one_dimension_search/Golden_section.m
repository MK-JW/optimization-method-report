clc
clear
%%  function test
alpha_lower = 0;
alpha_upper = 1;
tolerance = 10^-4;
x_current = -0.5;
d_current = 1;
[alpha_star] = golden_section(@f_test,alpha_lower,alpha_upper,tolerance,x_current,d_current);
%%  golden_section algorithm
function [alpha_star,x_next,f_next] = golden_section(f_test,alpha_lower,alpha_upper,tolerance,x_current,d_current)
    k = 0;
    alpha_lower_k = alpha_lower;
    alpha_upper_k = alpha_upper;
    reduction_rate = (sqrt(5) - 1)/2;
    l_k = alpha_upper_k - alpha_lower_k;
    alpha_left_k = alpha_upper_k - l_k*reduction_rate;
    alpha_right_k = alpha_lower_k + l_k*reduction_rate;
    %   k>=1时,开始迭代
    while((alpha_right_k - alpha_left_k)> tolerance)
        x_left = x_current + alpha_left_k*d_current;
        x_right = x_current + alpha_right_k*d_current;
        f_left = f_test(x_left);
        f_right = f_test(x_right);
        if(f_left> f_right)
            alpha_lower_k = alpha_left_k;
            l_k = alpha_upper_k - alpha_lower_k;
            alpha_left_k = alpha_upper_k - l_k*reduction_rate;
            alpha_right_k = alpha_lower_k + l_k*reduction_rate; 
        else
            alpha_upper_k = alpha_right_k;
            l_k = alpha_upper_k - alpha_lower_k;
            alpha_left_k = alpha_upper_k - l_k*reduction_rate;
            alpha_right_k = alpha_lower_k + l_k*reduction_rate;
        end
        k = k + 1;
    end
    alpha_star = (alpha_left_k+alpha_right_k)/2;
    x_next = x_current + alpha_star*d_current;
    f_next = f_test(x_next);
end
%%  function
function f_test = f_test(x)
f_test = 2*x^2 - x - 1;
end