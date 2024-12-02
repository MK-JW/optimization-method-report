clc;
clear;
format short
%%  给定初始参数
a_low = 0; a_up = 1; %给一个上下界
tol = 10^-8;
d_disturbance = 0.1*tol;
f = @(x1) 2*x1^2-x1-1;
k_count = 0;
%%  开始进行迭代
while(a_up - a_low>tol)
    a_middle = (a_up + a_low)/2;
    a1 = a_middle - d_disturbance;
    a2 = a_middle + d_disturbance;
    %   计算f(x+a1*d)与f(x+a2*d)
    if(f(a2)>f(a1))
        a_up = a2;
    elseif(f(a2)<f(a1))
        a_low = a1;
    else
        a_up = a2;
        a_low = a1;
    end
    k_count = k_count + 1;
end
%%  得到最后的步长
a = (a_up + a_low)/2;