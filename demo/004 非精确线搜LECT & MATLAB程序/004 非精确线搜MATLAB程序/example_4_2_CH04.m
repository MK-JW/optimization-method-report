function example_4_2_CH04
%-----------------------------------------------
%���˸�-Ӧ�����Ż�������MATLABʵ��CH04-4.2
%��4.2  Ŀ�꺯��f(x)=x1^2+x2^2-1;
%  x=[2 2]';  d=[-1 -1]'
%------------------------------------------------
x_current=[2;2];
d_current=[-1;-1];
rho=0.1;
[alpha_acceptable,x_next,f_next,k]=Armijo_search(@f_test2,@g_test2,x_current,d_current,rho)
%������f(x+alpha*d)����alpha��[0,9]����
N=100;
alpha_1=9;alpha_0=0;
da=(alpha_1-alpha_0)/N;
alpha_c=alpha_0:da:alpha_1;
for i=1:N+1
   xx(1,i)=x_current(1)+alpha_c(i)*d_current(1);
   xx(2,i)=x_current(2)+alpha_c(i)*d_current(2);
end
f=xx(1,:).^2+xx(2,:).^2-1;
plot(alpha_c,f,'k');
hold on
plot(alpha_acceptable,f_next,'*r');
%-------------------------------------------
%������f(x+alpha*d)����alpha��������alpha=0��
%����б��k=g(x)'*d
g_c=g_test2(x_current);
k=(g_c')*d_current;
f_c=f_test2(x_current);
NN=50;
ty_c=f_c+k*alpha_c(NN);
plot([alpha_c(1),alpha_c(NN)],[f_c,ty_c],'b');
%-----------------------------------------------
%��б��Ϊk_r=��g(x)'*d��*rho��alpha=0����ֱ��
k_r=k*rho;
NN=65;
ty_r=f_c+k_r*alpha_c(NN);
plot([alpha_c(1),alpha_c(NN)],[f_c,ty_r],'m');
end
function f_test2=f_test2(x)
x1=x(1);x2=x(2);
f_test2=x1^2+x2^2-1;
end
function g_test2=g_test2(x)
x1=x(1);x2=x(2);
g_test2=[2*x1;2*x2];
end