function example_3_9_CH03
x_current=3;
d_current=-1;
alpha_lower=0;
alpha_upper=2;
tolerance=1e-4;
[alpha_star,x_next,f_next,k]=Fibonacci_search(@f_test3,x_current,d_current,alpha_lower,alpha_upper,tolerance)
end

function f_test3=f_test3(x)
if(x<=2)
    f_test3=-x+3;
else
    f_test3=0.5*x;
end
end