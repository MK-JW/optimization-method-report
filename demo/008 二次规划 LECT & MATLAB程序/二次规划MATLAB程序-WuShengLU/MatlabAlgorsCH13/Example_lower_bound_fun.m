close;clear;clc;

figure(1)
t0=0;t1=pi*(1-1/6);
N=40;dt=(t1-t0)/N;
t=t0:dt:t1;
f=-5*sin(t);
plot(t,f)
hold on
x=[t0,pi/6,11*pi/20,4*pi/6,t(end)]';
fx=-5*sin(x);
plot(x,fx,'ob')
plot(x,-6*ones(5,1),'.b')
plot([-1,-1],[1,-6],'k')
plot([-1,4],[-6,-6],'k')
for i=1:5
    hold on
    plot([x(i),x(i)],[-6,fx(i)],'--b')
end
%axis padded
axis off
xx=x(2:4);
gi=-5*cos(xx);fi=-5*sin(xx);
n=length(t);
ffl=zeros(n,1);
for i=1:n
    fflk=fi+gi.*(t(i)*ones(length(xx),1)-xx);
    ffl(i)=max(fflk);
end
plot(t,ffl,'m')


mm=[3,4,6,10];
for ii=1:4
figure(2)
subplot(2,2,ii)
t0=pi/100;t1=pi*(1-1/6);
N=40;dt=(t1-t0)/N;
t=t0:dt:t1;
f=-5*sin(t);
plot(t,f,'b')
hold on
m=mm(ii);
dx=(t1-t0)/m;
x=(t0:dx:t1)';
fx=-5*sin(x);
plot(x,fx,'.b')
plot(x,-6*ones(length(x),1),'.b')
plot([-1,-1],[1,-6],'k')
plot([-1,4],[-6,-6],'k')
for i=1:length(x)
    hold on
    plot([x(i),x(i)],[-6,fx(i)],'-.k')
end
%axis padded
axis off
xx=x(2:(length(x)-1));
gi=-5*cos(xx);fi=-5*sin(xx);
n=length(t);
ffl=zeros(n,1);
for i=1:n
    fflk=fi+gi.*(t(i)*ones(length(xx),1)-xx);
    ffl(i)=max(fflk);
end
plot(t,ffl,'m')
str0=num2str(m-2);str=['k = ',str0]
title(str)
end

