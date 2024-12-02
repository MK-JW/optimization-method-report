function example_3_1_CH03
x_current=-0.5;
d_current=1;
alpha_lower=0;
alpha_upper=1;
tolerance=1e-4;
[alpha_star,x_next,f_next,k]=Dichotomous_search(@f_test1,x_current,d_current,alpha_lower,alpha_upper,tolerance)
end

function f_test1=f_test1(x)
f_test1=2*x^2-x-1;
end