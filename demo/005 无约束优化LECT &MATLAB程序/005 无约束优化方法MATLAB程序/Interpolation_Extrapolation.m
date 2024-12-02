close all;%interpolation and extrapolation
clear all;
clc;
f=@(t) t.*log(t);
g=@(t) 1+log(t);
N=100;
figure(1)
tL=0.005;
tU=1.2;
dt=(tU-tL)/N;
t=tL:dt:tU;
plot(t,f(t),'k')
hold on
% interpolation in [t0,0.9]
t1=tL;
t2=0.9;
plot(t1,f(t1),'or')
text(t1+0.003,f(t1)+0.3,'\alpha_1')
plot(t2,f(t2),'or')
text(t2+0.008,f(t1)+0.3,'\alpha_2')
t_star=t1-g(t1)*(t2-t1)^2/(2*(f(t2)-f(t1)-g(t1)*(t2-t1)))
plot(t_star,f(t_star),'*r');
text(t_star+0.02,f(t1)+0.35,'\alpha^*')
t_g=[t1,t_star,t2];
gt_t=g(t_g)
stem(t_g,gt_t,'ob')
title('\phi (\alpha) = \alpha*log(\alpha)')
xlabel('\alpha');ylabel('\phi (\alpha)   /   \phi^{,}(\alpha)');
axis ([-0.05,1.25,-5,1.2])

%extrapolation in [t0,0.1]
figure(2)
N=100;
tL=0.005;
tU=0.7;
dt=(tU-tL)/N;
t=tL:dt:tU;
plot(t,f(t),'k')
hold on
t1=tL;
t2=0.1;
plot(t1,f(t1),'or')
text(t1+0.003,f(t1)+0.3,'\alpha_1')
plot(t2,f(t2),'or')
text(t2,f(t1)+0.3,'\alpha_2')
t_star=t2-g(t2)*(t2-t1)/(g(t2)-g(t1))
plot(t_star,f(t_star),'*r');
text(t_star,f(t1)+0.35,'\alpha^*')
t_g=[t1,t_star,t2];
gt_t=g(t_g)
stem(t_g,gt_t,'ob')
title('\phi (\alpha) = \alpha*log(\alpha)')
xlabel('\alpha');ylabel('\phi (\alpha)   /   \phi^{,}(\alpha)');
axis ([-0.05,0.72,-5,1.2])

figure(3)
%subplot(2,1,1)
tL=0.005;
tU=1.2;
dt=(tU-tL)/N;
t=tL:dt:tU;
plot(t,f(t),'k')
hold on
% interpolation in [t0,0.9]
t1=tL;
t2=0.9;
plot(t1,f(t1),'or')
text(t1+0.003,f(t1)+0.3,'\alpha_1')
plot(t2,f(t2),'or')
text(t2+0.008,f(t1)+0.3,'\alpha_2')
t_star=t1-g(t1)*(t2-t1)^2/(2*(f(t2)-f(t1)-g(t1)*(t2-t1)))
plot(t_star,f(t_star),'*r');
text(t_star+0.02,f(t1)+0.35,'\alpha^*')
t_g=[t1,t_star,t2];
gt_t=g(t_g)
stem(t_g,gt_t,'ob')
title('\phi (\alpha) = \alpha*log(\alpha)')
xlabel('\alpha');ylabel('\phi (\alpha)   /   \phi^{,}(\alpha)');
axis ([-0.05,1.25,-5,1.2])

%subplot(2,1,2)
%tL=0.005;
%tU=0.7;
%dt=(tU-tL)/N;
%t=tL:dt:tU;
%plot(t,f(t),'k')
hold on
t1=tL;
t2=0.1;
%plot(t1,f(t1),'or')
text(t1+0.001,f(t1)+0.3,'\alpha_1')
plot(t2,f(t2),'og')
text(t2-0.004,f(t1)+0.3,'\alpha_2','Color','g')
t_star=t2-g(t2)*(t2-t1)/(g(t2)-g(t1))
plot(t_star+0.001,f(t_star),'*g');
text(t_star,f(t1)+0.35,'\alpha^*','Color','g')
t_g=[t_star,t2];
gt_t=g(t_g)
stem(t_g,gt_t,'og')
%title('\phi (\alpha) = \alpha*log(\alpha)')
%xlabel('\alpha');ylabel('\phi (\alpha)   /   \phi^{,}(\alpha)');
%axis ([-0.05,0.72,-5,1.2])
