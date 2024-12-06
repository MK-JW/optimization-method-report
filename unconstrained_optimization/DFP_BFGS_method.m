clc
clear
%%  测试写的function
%%  function DFP algorithm
function [x_optimal,f_optimal,k] = DFP_method(f_test,g_test,x_initial,tolerance)
k = 0;
x_current = x_initial;
n = length(x_current);
Q_current = eye(n);
f_current = f_test(x_current);
g_current = g_test(x_current);
k = 1;
x_next = x_current - Q_current*g_current;
f_next = f_test(x_next);
g_next = g_test(x_next);
while(abs(f_next - f_current)> tolerance)
    s_current = x_next - x_current;
    y_current = g_next - g_current;
    x_current = x_next;
    Q_next = Q_current + (s_current)*(s_current)'/(s_current)'*(y_current) - (Q_current*y_current)*(Q_current*y_current)'...
        /(y_current'*Q_current*y_current);  %DFP更新公式
    x_next = x_current - Q_next*g_next;    
    f_next = f_test(x_next);
    g_next = g_test(x_next);
    k = k+1;
end
end
%%  function BFGS algorithm
function [] = BFGS_method()
end
%%  function test