function example_3_5_CH03
x_current=[2,2];
d_current=[-1,-1];
alpha_lower=0;
alpha_upper=2;
tolerance=1e-6;
[alpha_star,x_next,f_next,k]=Trisection_search(@f_test2,x_current,d_current,alpha_lower,alpha_upper,tolerance)
end

function f_test2=f_test2(x)
x1=x(1);x2=x(2);
f_test2=x1^2+x2^2-1;
end