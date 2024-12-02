clc;
clear;
format short
%%  给定初始参数
a_low = -10; a_up = 10^3;
tol = 10^-3;
Fn = (a_up - a_low)/tol;
k_count = 0;
f = @(x) (1-x)^2;
%% Fibonacci数列求迭代数k_count
f0 = 1;
f1 = 1;
f2 = f0+f1;
while(Fn>f2)
    f0 = f1;
    f1 = f2;
    f2 = f0+f1;
    k_count = k_count+1;
end
%%  迭代计算精确步长
t0_reductionrate = f1/f2; %计算第一次迭代缩减率
for i = 1:k_count
    f2 = f1;
    f1 = f0;
    f0 = f2 - f1;
    t_reductionrate = f1/f2; %计算缩减率
    if(i == 1)
        a1 = a_up - (a_up - a_low)*t0_reductionrate;
        a2 = a_low + (a_up - a_low)*t0_reductionrate;
    end
    f_a1 = f(a1);
    f_a2 = f(a2);
    if(f_a1<f_a2)
        a_up = a2;
        a2 = a1;
        a1 = a_up - (a_up - a_low)*t_reductionrate;
    elseif(f_a1>f_a2)
        a_low = a1;
        a1 = a2;
        a2 = a_low + (a_up - a_low)*t_reductionrate;
    else
        a_low = a1;
        a_up = a2;
        a1 = a_up - (a_up - a_low)*t_reductionrate;
        a2 = a_low + (a_up - a_low)*t_reductionrate;
    end
end
% %%  迭代计算精确步长
% for i = 1:k_count
%     k = 1;
%     t_reductionrate = f1/f2; %计算缩减率
%     a1 = a_up - (a_up - a_low)*t_reductionrate;
%     a2 = a_low + (a_up - a_low)*t_reductionrate;
%     if(f(a1)<f(a2))
%         a_up = a2;
%     elseif(f(a1)>f(a2))
%         a_low = a1;
%     else
%         a_low = a1;
%         a_up = a2;
%     end
%     f2 = f1;
%     f1 = f0;
%     f0 = f2 - f1;
% end
%% 一次对分搜索
d_disturbance = 0.1*tol;
a_middle = (a_up+a_low)/2;
a1_new = a_middle - d_disturbance;
a2_new = a_middle + d_disturbance;
if(f(a1_new)<f(a2_new))
    a_up = a2_new;
elseif(f(a1_new)>f(a2_new))
    a_low = a1_new;
else
    a_low = a1_new;
    a_up = a2_new;
end
%%  得到最终结果
a = (a_low+a_up)/2;