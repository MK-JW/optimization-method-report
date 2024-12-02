function example_3_2_CH03
x_current=[2,2];
d_current=[-1,-1];
alpha_lower=0;
alpha_upper=2;
tolerance=1e-6;
%==============plot==================
N=100;N1=N+1;
d_a=(alpha_upper-alpha_lower)/N;
a=alpha_lower:d_a:alpha_upper;
f_a=ones(N1,1);
for i=1:N1
    x_a=x_current+a(i)*d_current;
    f_a(i,1)=f_test1(x_a);
end
plot(a,f_a(:,1),'k');
hold on
%===========================================   


[alpha_star,x_next,f_next,k]=Dichotomous_search(@f_test2,x_current,d_current,alpha_lower,alpha_upper,tolerance)
end

function f_test2=f_test2(x)
x1=x(1);x2=x(2);
f_test2=x1^2+x2^2-1;
end