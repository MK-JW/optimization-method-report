clear all;
clc;
ruo=0;
t1=0.3:0.05:0.7;
t2=0.3:0.05:0.7;
[x1,x2]=meshgrid(t1,t2);
f=(x1-1).^2+(x2-1).^2+ruo*(x1+x2-1).^2;
hold on
mesh(x1,x2,f,'EdgeColor','g');

hold on
ruo=50;
f=(x1-1).^2+(x2-1).^2+ruo*(x1+x2-1).^2;
mesh(x1,x2,f,'EdgeColor','b');

hold on
ruo=100;
f=(x1-1).^2+(x2-1).^2+ruo*(x1+x2-1).^2;
mesh(x1,x2,f,'EdgeColor','m');

x=0.5;y=0.5;
z=(x-1)^2+(y-1)^2;
plot3(x,y,z,'ko')