close;clear;clc;
figure(1)
hold on
tx=-3:2:3;ty=-3:2:3;
[x,y]=meshgrid(tx,ty);
z=1-y;
mesh(x,y,z)
plot3(-2,-2,3,'*r')
plot3(-2,-2,3,'or')
xlabel('x_1');ylabel('x_2');zlabel('x_3');
grid on


c=-5;%f_star=-18.8
f2=@(x,y,z) 0.5*x.^2+0.5*y.^2+2*x+y-z-c;
fimplicit3(f2,'EdgeColor','none','FaceAlpha',.4,'FaceColor',[0.3010 0.7450 0.9330])

hold on
c=-3;
f2=@(x,y,z) 0.5*x.^2+0.5*y.^2+2*x+y-z-c;
fimplicit3(f2,'EdgeColor','none','FaceAlpha',.4,'FaceColor',[0.4660 0.6740 0.1880])
c=-7;
f2=@(x,y,z) 0.5*x.^2+0.5*y.^2+2*x+y-z-c;
fimplicit3(f2,'EdgeColor','none','FaceAlpha',.4,'FaceColor',[0 0.4470 0.7410])




