function test_fmincon_equal
fun=@(x)x(1)^2+x(2)^2;
x0=rand(2,1);
A=[];   
Aeq=[];
b=[];
beq=[];
lb=[];
ub=[];
exitflag=1;
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
[x,fval,exitflag,output,lambda] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@mycon,options);
lambda.eqlin
lambda.eqnonlin
lambda.ineqlin
lambda.lower
lambda.upper
lambda.ineqnonlin
end

function [c,ceq]=mycon(x)
c=[];  %此处不要忘记将不等式改成不等式<=0的标准形式
ceq=x(1).*x(2)-1;
end