clc
clear
%%  一个例子
x_current1 = 5;d_current1 = 1;roh1 = 0.1;sigma1 = 0.1;
x_current2 = -2;d_current2 = 1;roh2 = 0.1;sigma2 = 0.1;
[alpha_acceptable1] = Armijo_wolfe_search(@f_test1,@g_test1,x_current1,d_current1,roh1,sigma1);
[alpha_acceptable2] = Armijo_wolfe_search(@f_test2,@g_test2,x_current2,d_current2,roh2,sigma2);
%%  Armijo_wolfe_search
function [alpha_acceptable] = Armijo_wolfe_search(f_test,g_test,x_current,d_current,rho,sigma)
    k_max = 1000;
    k = 0;
    alpha_lower_k = 0;
    alpha_upper_k = 10^8;
    x_alpha_lower_k = x_current + alpha_lower_k*d_current;
    f_x_alpha_lower_k = f_test(x_alpha_lower_k);
    df_x_alpha_lower_k = d_current'*g_test(x_alpha_lower_k);
    f_x_alpha_lower_0 = f_x_alpha_lower_k;
    df_x_alpha_lower_0 = f_x_alpha_lower_k;
    tolerance = 10^-15;
%     if(df_x_alpha_lower_0> 0)
%         df_x_alpha_lower_0 = -df_x_alpha_lower_0;
%     end
    if(abs(df_x_alpha_lower_k) >tolerance)
      alpha_k = -2*f_x_alpha_lower_k/df_x_alpha_lower_k;
    else
        alpha_k = 1;
    end
    if(alpha_k - alpha_lower_k <tolerance)
        alpha_k = 1;
    end
    for k = 1:k_max
        x_alpha_k = x_current + alpha_k*d_current;
        f_x_alpha_k = f_test(x_alpha_k);
        df_x_alpha_k = d_current'*g_test(x_alpha_k);
        Armijo_condition = f_x_alpha_k - f_x_alpha_lower_0 - rho*df_x_alpha_lower_0*alpha_k;
        wolfe_condition =  abs(df_x_alpha_k) - sigma*abs(df_x_alpha_lower_0);
        if(Armijo_condition <=0)
            if(wolfe_condition <=0)
                alpha_acceptable = alpha_k;
                break;
            else
                if(df_x_alpha_k <0)
                    delta_alpha_k = (alpha_k - alpha_lower_k)*df_x_alpha_k/(df_x_alpha_lower_k - df_x_alpha_k);
                    if(delta_alpha_k <=0)
                        alpha_k_temp = alpha_k - delta_alpha_k;
                    else
                        alpha_k_temp = alpha_k + delta_alpha_k;
                    end
                    alpha_lower_k = alpha_k;
                    f_x_alpha_lower_k = f_x_alpha_k;
                    df_x_alpha_lower_k = df_x_alpha_k;
                    alpha_k = alpha_k_temp;
                else
                    if(alpha_k<alpha_upper_k)
                        alpha_upper_k = alpha_k;
                    end
                    alpha_k_temp = alpha_lower_k - (1/2)*((alpha_k - alpha_lower_k)^2*df_x_alpha_lower_k)/(f_x_alpha_k - ...
                    f_x_alpha_lower_k - df_x_alpha_lower_k*(alpha_k - alpha_lower_k));
                    alpha_k = alpha_k_temp;
                end
            end
        else
            if(alpha_k <alpha_upper_k)
                alpha_upper_k = alpha_k;
            end
            alpha_k_temp = alpha_lower_k - (1/2)*((alpha_k - alpha_lower_k)^2*df_x_alpha_lower_k)/(f_x_alpha_k - ...
                f_x_alpha_lower_k - df_x_alpha_lower_k*(alpha_k - alpha_lower_k));
            alpha_k = alpha_k_temp;
        end
        if(alpha_upper_k - alpha_lower_k <tolerance)
            alpha_acceptable = alpha_k;
            break;
        end
    end
    if((Armijo_condition >0)||(wolfe_condition>0))
        alpha_acceptable = NaN;
    end
end
% %   二分法GPT写的
% function [alpha_acceptable] = Armijo_wolfe_search(f_test, g_test, x_current, d_current, rho, sigma)
%     k_max = 1000; % 最大迭代次数
%     tolerance = 1e-6; % 容差调整
%     alpha_lower_k = 0; % 步长下界
%     alpha_upper_k = 1e8; % 步长上界
%     alpha_k = 1; % 初始步长
% 
%     % 计算初始的函数值和梯度
%     f_x_alpha_lower_0 = f_test(x_current);
%     df_x_alpha_lower_0 = g_test(x_current) * d_current;
% 
%     % 迭代搜索
%     for k = 1:k_max
%         % 计算当前步长下的函数值和梯度
%         x_alpha_k = x_current + alpha_k * d_current;
%         f_x_alpha_k = f_test(x_alpha_k);
%         df_x_alpha_k = g_test(x_alpha_k) * d_current;
% 
%         % Armijo 条件：确保足够的下降
%         Armijo_condition = f_x_alpha_k - f_x_alpha_lower_0 - rho * df_x_alpha_lower_0 * alpha_k;
% 
%         % Wolfe 条件：确保梯度足够小
%         Wolfe_condition = abs(df_x_alpha_k) - sigma * abs(df_x_alpha_lower_0);
% 
%         % 如果同时满足 Armijo 和 Wolfe 条件，则接受当前步长
%         if (Armijo_condition <= 0 && Wolfe_condition <= 0)
%             alpha_acceptable = alpha_k;
%             disp(['Iteration ', num2str(k), ': alpha = ', num2str(alpha_k)]);
%             return;
%         end
% 
%         % 如果 Armijo 条件未满足（下降不足），则缩小步长
%         if (Armijo_condition > 0)
%             alpha_upper_k = alpha_k; % 缩小步长区间
%         else
%             alpha_lower_k = alpha_k; % 放大步长区间
%         end
% 
%         % 如果 Wolfe 条件未满足（梯度太大），则调整步长
%         if (Wolfe_condition > 0)
%             alpha_lower_k = alpha_k; % 放大步长区间
%         else
%             alpha_upper_k = alpha_k; % 缩小步长区间
%         end
% 
%         % 通过二分法更新步长
%         alpha_k = 0.5 * (alpha_lower_k + alpha_upper_k);
% 
%         % 如果步长区间收敛，退出
%         if (alpha_upper_k - alpha_lower_k < tolerance)
%             alpha_acceptable = alpha_k;
%             disp(['Iteration ', num2str(k), ': alpha = ', num2str(alpha_k)]);
%             return;
%         end
%     end
% 
%     % 如果超出最大迭代次数，则返回NaN
%     alpha_acceptable = NaN;
%     disp('Warning: Maximum iterations reached without finding acceptable alpha.');
% end
%%  例子函数1
function f_test1 = f_test1(x)
    f_test1 = (1/x)*sin(3*x);
end
%%  例子函数导数1
function g_test1 = g_test1(x)
    g_test1 = 3*cos(3*x)/x - sin(3*x)/x^2;
end
%%  例子函数2
function f_test2 = f_test2(x)
    f_test2 = -3*x*sin(0.75*x)+exp(-2*x);
end
%%  例子函数导数2
function g_test2 = g_test2(x)
    g_test2 = -2/exp(2*x)-3*sin((3*x)/4) - (9*x*cos((3*x)/4))/4;
end