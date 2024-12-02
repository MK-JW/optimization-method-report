function example_4_1_CH04
%-----------------------------------------------
%���˸�-Ӧ�����Ż�������MATLABʵ��CH04-4.2
%��4.1  Ŀ�꺯��f(x)=-s*x*sin(0.75*x)+exp(-2*x);
%  x=-2;  d=1
%------------------------------------------------
x_current=-2;
d_current=1;
rho=0.1;
[alpha_acceptable,x_next,f_next,k]=Armijo_search(@f_test1,@g_test1,x_current,d_current,rho);
%������f(x+alpha*d)����alpha��[0,9]����
N=100;
alpha_1=9;alpha_0=0;
da=(alpha_1-alpha_0)/N;
alpha_c=alpha_0:da:alpha_1;
xx=x_current+alpha_c*d_current;
f=-3*xx.*sin(0.75*xx)+exp(-2*xx);
plot(alpha_c,f,'k');
hold on
plot(alpha_acceptable,f_next,'*r');
%-------------------------------------------
%������f(x+alpha*d)����alpha��������alpha=0��
%����б��k=g(x)'*d
g_c=g_test1(x_current);
k=(g_c')*d_current;
f_c=f_test1(x_current);
NN=7;
ty_c=f_c+k*alpha_c(NN);
plot([alpha_c(1),alpha_c(NN)],[f_c,ty_c],'b');
%-----------------------------------------------
%��б��Ϊk_r=��g(x)'*d��*rho��alpha=0����ֱ��
k_r=k*rho;
NN=65;
ty_r=f_c+k_r*alpha_c(NN);
plot([alpha_c(1),alpha_c(NN)],[f_c,ty_r],'m');
end
function f_test1=f_test1(x)
f_test1=-3*x*sin(0.75*x)+exp(-2*x);
end
function g_test1=g_test1(x)
g_test1=-2/(exp(2*x))-3*sin((3*x)/4)-(9*x*cos((3*x)/4))/4;
end