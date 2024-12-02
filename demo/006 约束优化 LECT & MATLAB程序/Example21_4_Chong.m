function Example21_4_Chong
close;clear;clc;
fun=@(x) x(1).^2+x(2).^2+x(1).*x(2)-3*x(1);
x0=[1;1];
A=[]; b=[];  Aeq=[];beq=[];lb=[];ub=[];
exitflag=1;
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
[x,fval,exitflag,output,lambda]=fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@confun,options);

xx1=x(1);xx2=x(2);
f=@(x1,x2) x1.^2+x2.^2+x1.*x2-3*x1;
t=-1:0.1:6;  s=-4:0.1:4;
[X1,X2]=meshgrid(t,s);
F=f(X1,X2);
contour(X1,X2,F,'ShowText','on');
hold on
 
fx=f(xx1,xx2)
plot(xx1,xx2,'r*')
xxx1=xx1+0.2;xxx2=xx2-0.2;
fm=['f_m_i_n=',num2str(fx)];
text(xxx1,xxx2,fm)
plot([0,6],[0,0],'k')
plot([0,0],[0,4],'k')
xlabel('x_1');ylabel('x_2');
end
function [c,ceq]=confun(x)
c=[-x(1);-x(2)]; %inequality constraint g(x)<=0
ceq=[];% equality constraint
end

