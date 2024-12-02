clc;
clear;
format short
%%  给定一些初始参数
a_low = -10; a_up = 10^3;
tol = 10^-3;
k_count = 0;
f = @(x)(1-x)^2;
%%  进行迭代计算
while(a_up - a_low>tol)
    a_interval = (a_up - a_low)/4;
    a1 = a_low + a_interval;
    a2 = a_low + 2*a_interval;
    a3 = a_low + 3*a_interval;
    a_point = linspace(a_low,a_up,5);
    f_value = [f(a_low),f(a1),f(a2),f(a3),f(a_up)];
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
a = (a_low + a_up)/2;