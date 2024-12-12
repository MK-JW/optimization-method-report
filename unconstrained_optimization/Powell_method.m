clc
clear
addpath('D:\matlab\Optimization method\one_dimension_search');
%%  function test
x_initial = [0;0]; tolerance = 10^-15;
[a,b,c] = powell_method(@f_test,x_initial,tolerance);
%%  function powell algorithm
function [x_optimal,f_optimal,k] = powell_method(f_test,x_initial,tolerance)
k = 2;
x_current = x_initial;
d_1 = [1;0];
alpha_star = Fibonacci(f_test,-5,5,10^-9,x_current,d_1);
x_next1 = x_current +  alpha_star*d_1;
disp(x_next1);
d_2 = [0;1];
alpha_star = Fibonacci(f_test,-5,5,10^-9,x_next1,d_2);
x_next2 = x_next1 + alpha_star*d_2;
disp(x_next2);
while(norm(x_next2 - x_current)> tolerance)
    k = k+2;
    x = 2*x_next2 - x_current;
    if(f_test(x - x_next1)> f_test(x - x_current))
        d_3 = x - x_current;
    else
        d_3 = x - x_next1;
    end
    x_current = x_next2;
    d_1 = d_2;
    d_2 = d_3;
    alpha_star = Fibonacci(f_test,-5,5,10^-9,x_current,d_1);
    x_next1 = x_current +  alpha_star*d_1;
    alpha_star = Fibonacci(f_test,-5,5,10^-9,x_next1,d_2);
    x_next2 = x_next1 + alpha_star*d_2;
end
x_optimal = x_next2;
f_optimal = f_test(x_optimal);
end
%% test function 
function f_test = f_test(x)
x1 = x(1);
x2 = x(2);
f_test = (x1-2)^2 + (x2+3)^2;
end