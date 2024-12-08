clc
clear
addpath('D:\matlab\Optimization method\one_dimension_search');
%%  测试写的function
x_initial = [1;1]; tolerance = 10^-15;
[a,b,c] = DFP_method(@f_test1,@g_test1,x_initial,tolerance);
%%  function DFP algorithm
function [x_optimal,f_optimal,k] = DFP_method(f_test,g_test,x_initial,tolerance)
% % %   这里面的步长和牛顿法一样
% k = 1;
% x_current = x_initial;
% n = length(x_current);
% Q_current = eye(n);
% f_current = f_test(x_current);
% g_current = g_test(x_current);
% x_next = x_current - Q_current*g_current;
% f_next = f_test(x_next);
% g_next = g_test(x_next);
% while(norm(x_next - x_current)> tolerance) %这里面计算的时候要取模 
%     k = k+1;
%     s_current = x_next - x_current;
%     y_current = g_next - g_current;
%     x_previous = x_current;
%     x_current = x_next;
%     f_current = f_test(x_current);
%     g_current = g_test(x_current);
%     Q_previous = Q_current;
%     Q_current = Q_previous + ((s_current)*(s_current)')/((s_current)'*(y_current)) - ((Q_previous*y_current)*(Q_previous*y_current)')...
%         /(y_current'*Q_previous*y_current);  %DFP更新公式
%     x_next = x_current - Q_current*g_current;    
%     f_next = f_test(x_next);
%     g_next = g_test(x_next);
%     x_optimal = x_next;
%     f_optimal = f_test(x_optimal);
% end
% end
%%  结合一下一维搜索方法
k = 1;
sigma = 0.11;
rho = 0.1;
x_current = x_initial;
n = length(x_current);
Q_current = eye(n);
g_current = g_test(x_current);
d_current = -Q_current*g_current;
[alpha_acceptable] = Fibonacci(f_test,0,2,10^-9,x_current,d_current);   %斐波那契额精确搜索
% [alpha_acceptable] = Armijo_wolfe_search(f_test,g_test,x_current,d_current,rho,sigma);% 强wolfe非精确搜索
x_next = x_current + alpha_acceptable*d_current;
f_next = f_test(x_next);
while(norm(x_next - x_current)> tolerance)
    k = k+1;
    x_previous = x_current;
    x_current = x_next;
    g_previous = g_test(x_previous);
    g_current = g_test(x_current);
    Q_previous = Q_current;
    s_current = x_current - x_previous;
    y_current = g_current - g_previous;
    if(s_current'*y_current<= 0)    %拟牛顿方向正定条件
        Q_current = eye(n);
    else
        Q_current = Q_previous + ((s_current)*(s_current)')/((s_current)'*(y_current)) - ((Q_current*y_current)*(Q_current*y_current)')...
        /(y_current'*Q_current*y_current);  %DFP更新公式
    end
    d_current = -Q_current*g_current;
%     [alpha_acceptable] = Armijo_wolfe_search(f_test,g_test,x_current,d_current,rho,sigma);
    [alpha_acceptable] = Fibonacci(f_test,0,2 ,10^-9,x_current,d_current);
    if(isnan(alpha_acceptable)) %放宽Armijo_wolfe条件
        rho = rho + 0.1;
        sigma = sigma + 0.1;
        continue;
    else
        x_next = x_current + alpha_acceptable*d_current;
        f_next = f_test(x_next);
    end
end
x_optimal = x_next;
f_optimal = f_next;
end
% %%  function BFGS algorithm
% function [] = BFGS_method()
% end
%%  function test
function f_test1 = f_test1(x)
    x1 = x(1);
    x2 = x(2);
    f_test1 =2* x1^2 + x2^2 - x1 + 3;
end
function g_test1 = g_test1(x)
    x1 = x(1);
    x2 = x(2);
    g1 = 4*x1 - 1;
    g2 = 2*x2;
    g_test1 = [g1;g2];
end