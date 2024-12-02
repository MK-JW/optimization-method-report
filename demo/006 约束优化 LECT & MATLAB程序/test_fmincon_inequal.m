function test_fmincon_inequal
fun=@(x)(x(1)-1).^2+(x(2)-2).^2;
%x0=rand(2,1);
x0=[-3;0];
A=[];   
Aeq=[];
b=[];
beq=[];
lb=[];
ub=[];
exitflag=1;
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
[x,fval,exitflag,output,lambda]=fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@mycon,options);
x
lambda.eqlin
lambda.eqnonlin
lambda.ineqlin
lambda.lower
lambda.upper
lambda.ineqnonlin

s=-4:0.02:4;t=-1:0.02:5;
[u,v]=meshgrid(s,t);
f=(u-1).^2+(v-2).^2;
figure(1)
contour(u,v,f,'LineWidth',0.9)
hold on
plot(1,2,'ok')
g1=-u+1-(v-2).^2;
vc=[0,0];
contour(u,v,g1,vc,'Color','k','LineWidth',0.9)
g2=-u+1-(v-2).^2+1;
contour(u,v,g2,vc,'Color','g','LineWidth',0.9)
g3=-u+1-(v-2).^2-1;
contour(u,v,g3,vc,'Color','r','LineWidth',0.9)


end

function [c,ceq]=mycon(x)
% c=-(x(1)-1)-(x(2)-2).^2;  %inequality constraint g(x)<=0
% c=-(x(1)-1)-(x(2)-2).^2+1; 
  c=-(x(1)-1)-(x(2)-2).^2-1; 
ceq=[];% equality constraint
end
