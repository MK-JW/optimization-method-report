clc
clear
addpath('D:\matlab\Optimization method\one_dimension_search');
%%  test functon
x_initial = [1;1]; tolerance = 10^-15;
[a,b,c] = negative_gradient(@f_test,@g_test,x_initial,tolerance);
%%  Negative_gradient_method
function [x_optimal,f_optimal,k] = negative_gradient(f_test,g_test,x_initial,tolerance)
k = 1;
rho = 0.1;
sigma = 0.11;
x_current = x_initial;
g_current = g_test(x_current);
d_current = -g_current;
[alpha_acceptable] = Armijo_wolfe_search(f_test,g_test,x_current,d_current,rho,sigma);  %强wolfe条件
% [alpha_star] = Fabonacci(f_test,0,2,10^-9,x_current,d_current); %斐波那契
x_next = x_current + alpha_acceptable*d_current;
% x_next = x_current + alpha_star*d_current;
while(norm(x_next - x_current)> tolerance)
    k = k+1;
    x_current = x_next;
    g_current = g_test(x_current);
    d_current = -g_current;
    x_next = x_current + alpha_acceptable*d_current;
%     x_next = x_current + alpha_star*d_current;
end
x_optimal = x_next;
f_optimal = f_test(x_optimal);
end
%%  function
function f_test = f_test(x)
x1 = x(1);
x2 = x(2);
f_test = 2*x1^2 + x2^2 - x1 + 3;
end

function g_test = g_test(x)
x1 = x(1);
x2 = x(2);
g1 = 4*x1 - 1;
g2 = 2*x2;
g_test = [g1;g2];
end