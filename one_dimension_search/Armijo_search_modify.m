clc
clear
%%  测试用的函数
x_current = 9;
d_current = 1;
roh = 10^-1;
[alpha_acceptable] = Armijo_search(@f_test1,@g_test1,x_current,d_current,roh);
%%  Armjo条件搜索非精确步长
function [alpha_acceptable]=Armijo_search(f_test,g_test,x_current,d_current,rho)
k_max = 1000;
k = 0;
alpha_lower_k = 0;
alpha_upper_k = 10^8;
x_alpha_lower_k = x_current + alpha_lower_k*d_current;
f_x_alpha_lower_k = f_test(x_alpha_lower_k);
df_x_alpha_lower_k = g_test(x_alpha_lower_k)*d_current;     %链式法则得到关于alpha的斜率
% alpha_k = alpha_upper_k;
alpha_k = -f_x_alpha_lower_k/g_test(x_alpha_lower_k);
for k = 1:k_max
    x_alpha_k = x_current + alpha_k*d_current;
    f_x_alpha_k = f_test(x_alpha_k);
    Armijo_condition = f_x_alpha_k - df_x_alpha_lower_k*rho*alpha_k - f_x_alpha_lower_k;
    if(Armijo_condition <=0)
        alpha_acceptable = alpha_k;
        break;
    else
        if(alpha_k <alpha_upper_k)
            alpha_upper_k = alpha_k;
        end
        alpha_k = alpha_lower_k - (1/2)*((alpha_k - alpha_lower_k)^2*df_x_alpha_lower_k)/(f_x_alpha_k - f_x_alpha_lower_k - ...
            df_x_alpha_lower_k*(alpha_k - alpha_lower_k));
%         alpha_k = alpha_lower_k + (1/2)*((alpha_k - alpha_lower_k)^2*df_x_alpha_lower_k)/(f_x_alpha_lower_k ...
%             - f_x_alpha_k + df_x_alpha_lower_k*(alpha_k - alpha_lower_k));
    end
end
if(k == k_max)
    alpha_acceptable = NaN;
end
end
%%  给定一个函数
function f_test1 = f_test1(x)
    f_test1 = (1/x)*sin(3*x);
end
%%  给定函数的梯度
function g_test1 = g_test1(x)
    g_test1 = 3*cos(3*x)/x - sin(3*x)/x^2;
end