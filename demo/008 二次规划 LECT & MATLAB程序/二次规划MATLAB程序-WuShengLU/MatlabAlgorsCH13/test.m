function test
x0=[0;0];
A0=[16,0;0,9];
plot(x0(1),x0(2),'r+')
hold on
h1=ezplot(@(x,y)ellipse(x,y,x0,A0),[-8,12,-8,12]);
set(h1,'color','b','linestyle',':'); % E_k-1
end

function z = ellipse(x,y,x0,A)
t=A(1,1)*A(2,2)-A(2,1)*A(1,2);
a=A(2,2)/t;d=A(1,1)/t;b=-A(1,2)/t;c=-A(2,1)/t;
z=a*(x-x0(1)).^2+(b+c)*(x-x0(1)).*(y-x0(2))+d*(y-x0(2)).^2-1;
end