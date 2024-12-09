clc
clear;
addpath('D:\matlab\Optimization method\one_dimension_search');
%%  测试函数
x_initial = [1;1]; tolerance = 10^-12;
[m,j,w] = conjugate_gradient(@f_test,@g_test,x_initial,tolerance);
%%  Conjugate_gradient algorithm
function  [x_optimal,f_optimal,k] = conjugate_gradient(f_test,g_test,x_initial,tolerance)
k  = 1;
rho = 0.1; sigma = 0.11;
x_current = x_initial;
n = length(x_current);
g_current = g_test(x_current);
d_current = -g_current;
[alpha_acceptable] = Armijo_wolfe_search(f_test,g_test,x_current,d_current,rho,sigma);
x_next = x_current + alpha_acceptable*d_current;
while(norm(x_next - x_current)> tolerance)
    k  = k+1;
    g_previous = g_current;
    d_previous = d_current;
    x_current = x_next;
    g_current = g_test(x_current);
    if(mod(k,n) == 0)
        d_current = -g_current;
    else
        beta_current = (g_current'*g_current)/(d_previous'*(g_current - g_previous));
        d_current = -g_current + beta_current*d_previous;
    end
    [alpha_acceptable] = Armijo_wolfe_search(f_test,g_test,x_current,d_current,rho,sigma);
    x_next = x_current + alpha_acceptable*d_current;
end
x_optimal = x_next;
f_optimal = f_test(x_optimal);
end
%%  function test
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