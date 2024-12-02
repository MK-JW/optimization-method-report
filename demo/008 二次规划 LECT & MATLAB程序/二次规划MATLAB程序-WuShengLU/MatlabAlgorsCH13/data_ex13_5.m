function data_ex13_5
% data_ex13_5
%initial point X0=[x0,L0]'
x0=[1.5,0.5,2.5,4.0]';
L0=1;
X0=[x0;L0];
epsi = 1e-7;

x=X0
c = cex_13_5(x)
h = hex_13_5(x)
end

%================cex_13_5.m============================================
% Program: cex_13_5.m
% Description: This function can be used to 
% test function kelley_ie.
% x=[x_points;L],example: x=X0=[x0;L0]
 function c = cex_13_5(x)
 nx=length(x);
 L = x(nx);
 xx = x(1:4);
 H=[1,0,-2,0;0,1,0,-2;-2,0,1,0;0,-2,0,1];
 p=zeros(nx-1,1);
 Q1=-[0.5,0,0,0;0,2,0,0;0,0,0,0;0,0,0,0];
 Q2=-[0,0,0,0;0,0,0,0;0,0,5/4,3/4;0,0,3/4,5/4];
 q1=[0.5;0;0;0];q2=[0;0;11/2;13/2];
 r1=3/4;r2=-35/2;
 c1 = L - (0.5*xx'*H*xx +xx'*p);
 c2 = 0.5*xx'*Q1*xx +xx'*q1+r1;
 c3 = 0.5*xx'*Q2*xx +xx'*q2+r2;
 c = [c1 c2 c3]';
 end

 %================hex_13_5.m============================================
 % Program: hex_13_5.m
% Description: This function can be used to 
% test function kelley_ie.
 function h = hex_13_5(x)
 nx=length(x);
 L = x(nx);
 xx = x(1:4);
 H=[1,0,-2,0;0,1,0,-2;-2,0,1,0;0,-2,0,1];
 p=zeros(nx-1,1);
 Q1=-[0.5,0,0,0;0,2,0,0;0,0,0,0;0,0,0,0];
 Q2=-[0,0,0,0;0,0,0,0;0,0,5/4,3/4;0,0,3/4,5/4];
 q1=[0.5;0;0;0];q2=[0;0;11/2;13/2];
 r1=3/4;r2=-35/2;

 h1 = [H*xx+p; -1];
 h2 = [-Q1*xx-q1;0];
 h3 = [-Q2*xx-q2;0];
 
 h = [h1 h2 h3];
 end