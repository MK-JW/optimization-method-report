clc;
clear;
format short
%%  给定一些初始参数
a_low = 0; a_up = 1;
tol = 10^-4;
k_count = 0;
x_current = -0.5;
d_current = 1;
f = @(x) 2*x^2 - x -1;
%%  进行迭代计算
while(a_up - a_low>tol)
    a_interval = (a_up - a_low)/4;
    a1 = a_low + a_interval;
    a2 = a_low + 2*a_interval;
    a3 = a_low + 3*a_interval;
    a_point = linspace(a_low,a_up,5);
    x_alow = x_current + a_low*d_current;
    x_a1 = x_current + a1*d_current;
    x_a2 = x_current + a2*d_current;
    x_a3 = x_current + a3*d_current;
    x_aup = x_current + a_up*d_current;
    f_value = [f(x_alow),f(x_a1),f(x_a2),f(x_a3),f(x_aup)];
    [~,I] = min(f_value);
    if(I == 1 )
        a_up = a_point(I+1);
    elseif(I == 5)
        a_low = a_point(I-1);
    else
        a_low = a_point(I-1);
        a_up = a_point(I+1);
    end
    k_count = k_count + 1;
end
%%  得到结果
alpha_star = (a_low + a_up)/2;