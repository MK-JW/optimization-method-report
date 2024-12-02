

[X,Y,Z] = peaks(40);
figure
%surfc(X,Y,Z)
meshc(X,Y,Z)
hold on
x0=-2.5;y0=3;
dx=1;dy=-2.8;
t=0:0.01:2.3;
N=length(t);
 x=x0+t*dx;y=y0+t*dy;
z =  3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ... 
   - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) ... 
   - 1/3*exp(-(x+1).^2 - y.^2);
plot3(x,y,z,'r');
xlabel('x1');ylabel('x2')
plot3(x(1),y(1),z(1),'*b');
z0=-10*ones(N,1);
plot3(x,y,z0,'r')
plot3(x(1),y(1),z0(1),'*b');