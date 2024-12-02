function example_3_2_3D_CH03()
close all;
clear all;
clc;

t1=-3:0.25:3;t2=-3:0.25:3;
[X,Y]=meshgrid(t1,t2);
Z=X.^2+Y.^2-1;
figure
%surfc(X,Y,Z)
meshc(X,Y,Z)
hold on
x0=2;y0=2;
dx=-1;dy=-1;
t=0:0.1:5;
N=length(t);
 x=x0+t*dx;y=y0+t*dy;
z = x.^2+y.^2-1;
plot3(x,y,z,'r');
xlabel('x1');ylabel('x2')
plot3(x(1),y(1),z(1),'*b');
z0=-5*ones(N,1);
plot3(x,y,z0,'r')
plot3(x(1),y(1),z0(1),'*b');

x_current=[2,2];
d_current=[-1,-1];
alpha_lower=0;
alpha_upper=2;
tolerance=1e-6;
[alpha_star,x_next,f_next,alpha_l,alpha_r,k]=Dichotomous_search1(@f_test2,x_current,d_current,alpha_lower,alpha_upper,tolerance);
x_l=x0+alpha_l*dx;
y_l=y0+alpha_l*dy;
z_l=x_l.^2+y_l.^2-1;
nzl=length(z_l);
plot3(x_l,y_l,z_l,'ob');
zl0=-5*ones(nzl,1);
plot3(x_l,y_l,zl0,'ob');
x_r=x0+alpha_r*dx;
y_r=y0+alpha_r*dy;
z_r=x_r.^2+y_r.^2-1;
nzr=length(z_r);
plot3(x_r,y_r,z_r,'pg');
zr0=-5*ones(nzr,1);
plot3(x_r,y_r,zr0,'pg');
plot3(x_next(1),x_next(2),f_next,'*r')
plot3(x_next(1),x_next(2),-5,'*r')
end

function f_test2=f_test2(x)
x1=x(1);x2=x(2);
f_test2=x1^2+x2^2-1;
end
function [alpha_star,x_next,f_next,alpha_l,alpha_r,k]=Dichotomous_search1(f_test,x_current,d_current,alpha_lower,alpha_upper,tolerance)

if(tolerance>=1e-8)
    disturbance_quantity=1e-9;
else
    disturbance_quantity=0.1*tolerance;
end

k=1;
alpha_lower_k=alpha_lower;alpha_l(k)=alpha_lower_k;
alpha_upper_k=alpha_upper;alpha_r(k)=alpha_upper_k;
%--------------------------------------------------------------------------
% k>=1时
%--------------------------------------------------------------------------
while (abs(alpha_upper_k-alpha_lower_k)>tolerance)
    k=k+1;
    alpha_middle_k=0.5*(alpha_upper_k+alpha_lower_k);
    alpha_left_k=alpha_middle_k-disturbance_quantity;alpha_l(k)=alpha_lower_k;
    alpha_right_k=alpha_middle_k+disturbance_quantity;alpha_r(k)=alpha_upper_k;
    x_alpha_left_k=x_current+alpha_left_k*d_current;
    x_alpha_right_k=x_current+alpha_right_k*d_current;
    f_alpha_left_k=f_test(x_alpha_left_k);
    f_alpha_right_k=f_test(x_alpha_right_k);

    if(f_alpha_left_k<f_alpha_right_k)
        alpha_upper_k=alpha_right_k;
    elseif(f_alpha_left_k>f_alpha_right_k)
        alpha_lower_k=alpha_left_k;
    else
        alpha_lower_k=alpha_left_k;
        alpha_upper_k=alpha_right_k;
    end
end
%--------------------------------------------------------------------------
% 取区间[alpha_lower_k，alpha_upper_k]的中点作为最佳步长
%--------------------------------------------------------------------------
alpha_star=0.5*(alpha_lower_k+alpha_upper_k);
x_next=x_current+alpha_star*d_current;
f_next=f_test(x_next);
end
