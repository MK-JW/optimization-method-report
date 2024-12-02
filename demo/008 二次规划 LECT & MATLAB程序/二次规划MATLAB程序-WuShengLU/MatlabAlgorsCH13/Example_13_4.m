function Example_13_4
%=============================================
% Example:  x0 mu0 nofeasible
% Find a minimizer of the QP problem
% minimize 0.5*x'*H*x +x'*p
% subject to A*x = b, x >= 0
% Solution:
% Execute the following commands:
% [H,p,A,b,x0,lmd0,mu0,rho,epsi] = data_ex13_3_2
% [xs,fs,k]= qp_path_nf(H,p,A,b,x0,lmd0,mu0,rho,epsi)
%==============================================
close;clear;clc;
%=====for plot3 for xk=================
global x_rec
%=======================================
[H,p,A,b,x0,lmd0,mu0,rho,epsi] = data_ex13_4;
[xs,fs,k]= qp_path_nf(H,p,A,b,x0,lmd0,mu0,rho,epsi)

xi=ones(4,1);xf=xi;
for i=1:4
    xi(i)=x0(i)-x0(4+i);
    xf(i)=xs(i)-xs(4+i);
end
r=[xf(1);xf(2)];s=[xf(3);xf(4)];
format short
r
s
dist=norm(r-s,2)
str=['the shortest distance between R and S = ',num2str(dist)]
title(str)

figure(1)
hold on
plot([0,4],[0,0],'k')
plot([0,0],[0,4],'k')
xlabel('x_1 & x_3');ylabel('x_2 & x_4');
plot([0,2,0,0],[0,0,1,0],'b');%R
plot([2,0,1,2],[2,3,2,2],'g');%S
tx=0.5;ty=0.25;
text(tx,ty,'R','Color','b')
tx=0.8;ty=2.35;
text(tx,ty,'S','Color','g')
plot(xi(1),xi(2),'or');plot(xi(3),xi(4),'om');
plot(xf(1),xf(2),'*r');plot(xf(3),xf(4),'*m');
plot([xf(1),xf(3)],[xf(2),xf(4)],'-.k')

[n0,m0]=size(x_rec);
xx=zeros(4,m0);
for j=1:m0
    for i=1:4
        xx(i,j)=x_rec(i,j)-x_rec(i+4,j);
    end
end
hold on
plot(xx(1,:),xx(2,:),'--r')
plot(xx(1,:),xx(2,:),'.r')
plot(xx(3,:),xx(4,:),'--g')
plot(xx(3,:),xx(4,:),'.m')
axis padded

kl=0:1:m0-1;
figure(2)
subplot(2,2,1)
plot(kl,xx(1,:),'b')
hold on
plot(kl,xx(1,:),'.r')
xlabel('Num of iter.');ylabel('x_1');
subplot(2,2,2)
plot(kl,xx(2,:),'b')
hold on
plot(kl,xx(2,:),'.r')
xlabel('Num of iter.');ylabel('x_2');
subplot(2,2,3)
plot(kl,xx(3,:),'b')
hold on
plot(kl,xx(3,:),'.r')
xlabel('Num of iter.');ylabel('x_3');
subplot(2,2,4)
plot(kl,xx(4,:),'b')
hold on
plot(kl,xx(4,:),'.r')
xlabel('Num of iter.');ylabel('x_4');

distt=(xx(1,:)-xx(3,:)).^2+(xx(2,:)-xx(4,:)).^2;
figure(3)
plot(kl,distt,'k')
hold on
plot(kl,distt,'ob')
xlabel('Num of iter.');ylabel('dist');


end
%
%
%
%===============================================================
function [H,p,A,b,x0,lmd0,mu0,rho,epsi] =data_ex13_4
%=============================================
% Example: 13.4   x0 is nofeasible point
% Find a minimizer of the QP problem
% minimize 0.5*x'*H*x +x'*p
% subject to A*x = b, x >= 0
%==============================================
% where
A0=[1,0,0,0;0,1,0,0;-1,-2,0,0;0,0,0,1;0,0,1,1;0,0,-1,-2];
n0=4;m=6;n=14;
p0=zeros(n0,1);
b=[0,0,-2,2,3,-6]';
H0=[1,0,-1,0;0,1,0,-1;-1,0,1,0;0,-1,0,1];
H00=zeros(n0,m);
H=[H0,-H0,H00;-H0,H0,H00;H00',H00',zeros(m,m)];
p=[p0;-p0;zeros(m,1)];
A=[A0,-A0,-eye(m,m)];
x0=ones(n,1);
lmd0=-ones(m,1);
mu0=ones(n,1);

format long
rho=n+20*sqrt(n);
epsi=1.e-5;
D=eig(H);
[n,m]=size(H);n=min(n,m);
if min(D)<=0
    H=H+1.e-9*eye(n,n);
end
end


%================================================================
% Program: qp_path_nf.m
% Title: Nonfeasible-initialization primal-dual 
% path-following algorithm for convex QP problems
% Description: Implements Algorithm 13.3 for convex QP 
% problems.
% Theory: See Practical Optimization Sec. 13.4.3.
% Input:     
%             H -- positive semidefinite Hessian matrix
%             A -- full row-rank constraint matrix A
%          p, b -- input vectors
% (x0,lmd0,mu0) -- initial point with x0 > 0 and mu0 > 0
%           rho -- parameter involved in determining 
%                  the value of tau, rho is required 
%                  to be no less than the square root of n
%          epsi -- tolerance for duality gap
% Output:   
%            xs -- solution vector
%            fs -- value of objective function at xs
%             k -- number of iterations at convergence
% Example:
% Find a minimizer of the QP problem
% minimize 0.5*x'*H*x +x'*p
% subject to A*x = b, x >= 0
% where
% H = [4 0 0; 0 1 -1; 0 -1 1]
% p = [-8 -6 -6]'
% A = [1 1 1]
% b = 3
% using the initial value
% x0 = [1 2 2]'
% Solution:
% Execute the following commands:
% H = [4 0 0; 0 1 -1; 0 -1 1]
% p = [-8 -6 -6]'
% A = [1 1 1]
% b = 3
% x0 = [1 2 2]'
% lmd0 = -1
% mu0 = [0.2 0.2 0.2]'
% rho = 3
% epsi = 1e-6
% [xs,fs,k]= qp_path_nf(H,p,A,b,x0,lmd0,mu0,rho,epsi)
% ===================================================
function [xs,fs,k]= qp_path_nf(H,p,A,b,x0,lmd0,mu0,rho,epsi)
%=====for plot3 for xk=================
global x_rec
%=======================================
disp(' ')
disp('Program qp_path_nf.m')
% Data initialization
n = length(x0);
rho_n = n + rho;
a_max = 1-1e-6;
x = x0(:);
%===========================================
x_rec=zeros(length(x0),50000);
k_rec=1;x_rec(:,k_rec)=x0;
%===========================================
lmd = lmd0(:);
mu = mu0(:);
xm = x.*mu;
gap = sum(xm);
k = 0;
% Iteration begins.
while gap > epsi
  % Compute Y0, yd, etc.
  tau = gap/rho_n;
  r_p = b - A*x;
  r_d = H*x + p - A'*lmd-mu;
  M = diag(mu);
  X = diag(x);
  Gam = inv(M+X*H);
  GaX = Gam*X*A';
  Y0 = inv(A*GaX);
  yd = Gam*(x.*(mu+r_d) - tau);
  % Calculate d_x,d_lmd, and d_mu.
  d_lmd = Y0*(A*yd + r_p);
  d_x = GaX*d_lmd - yd;
  d_mu = H*d_x - A'*d_lmd + r_d;
  % Calculate step size.
  ind = find(d_x < 0);
  a_p = min(x(ind)./(-d_x(ind)));
  ind = find(d_mu < 0);
  a_d = min(mu(ind)./(-d_mu(ind)));
  a_k = a_max*min([a_p a_d]);
  % Form new iterate. 
  x = x + a_k*d_x;
%===========================================
k_rec=k_rec+1;x_rec(:,k_rec)=x;
%===========================================
  lmd = lmd + a_k*d_lmd;
  mu = mu + a_k*d_mu;
  % Compute duality gap and check convergence.
  xm = x.*mu;
  gap = sum(xm);
  % Update iteration index, etc.
  k = k + 1;
end
xs = x;
fs = 0.5*xs'*(H*xs + 2*p);
%===========================================
x_rec=x_rec(:,1:k_rec);
%===========================================
end