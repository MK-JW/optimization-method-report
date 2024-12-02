function Example_13_3_2
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
[H,p,A,b,x0,lmd0,mu0,rho,epsi] = data_ex13_3_2;
[xs,fs,k]= qp_path_nf(H,p,A,b,x0,lmd0,mu0,rho,epsi)

figure(1)
hold on
xlabel('x_1');ylabel('x_2');zlabel('x_3');
grid on
%ellipsoid
c=-18.5;%f_star=-18.5
f=@(x,y,z) 2*x.^2+0.5*y.^2+0.5*z.^2-y.*z-8*x-6*y-6*z-c;
%interval = [-0.5 1.5 -0.5 1.5 -1 3];
%interval = [-1 3 -1 3 -1 3];
%fimplicit3(f,interval,'EdgeColor','none','FaceAlpha',.4,'FaceColor',[0.3010 0.7450 0.9330])
fimplicit3(f,'EdgeColor','none','FaceAlpha',.4,'FaceColor',[0.3010 0.7450 0.9330])
%x,y,z-axis
t=1.5;
plot3([0,t],[0,0],[0,0],'--k')
plot3([0,0],[0,t],[0,0],'--k')
plot3([0,0],[0,0],[0,t],'--k')
% Ax=b
tx=-0.5:0.5:1.5;ty=-0.5:0.5:1.5;
[x,y]=meshgrid(tx,ty);
z=3-x-y;
mesh(x,y,z)
% iterate points
plot3(x_rec(1,:),x_rec(2,:),x_rec(3,:),'-m')
plot3(x_rec(1,:),x_rec(2,:),x_rec(3,:),'.m')
plot3(x_rec(1,end),x_rec(2,end),x_rec(3,end),'*r')


x=x_rec(1,:);y=x_rec(2,:);z=x_rec(3,:);
ff=f(x,y,z)+c;
[n,m]=size(x_rec);
kl=0:1:m-1;
figure(2)
subplot(2,2,1)
plot(kl,x_rec(1,:),'b')
hold on
plot(kl,x_rec(1,:),'.r')
xlabel('Num of iter.');ylabel('x_1');
subplot(2,2,2)
plot(kl,x_rec(2,:),'b')
hold on
plot(kl,x_rec(2,:),'.r')
xlabel('Num of iter.');ylabel('x_2');
subplot(2,2,3)
plot(kl,x_rec(3,:),'b')
hold on
plot(kl,x_rec(3,:),'.r')
xlabel('Num of iter.');ylabel('x_3');
subplot(2,2,4)
kl=0:1:length(ff)-1;
plot(kl,ff,'k')
hold on
plot(kl,ff,'ob')
xlabel('Num of iter.');ylabel('f');

end
%

function [H,p,A,b,x0,lmd0,mu0,rho,epsi] = data_ex13_3_2
%=============================================
% Example: 13.3.2   x0 is nofeasible point
% Find a minimizer of the QP problem
% minimize 0.5*x'*H*x +x'*p
% subject to A*x = b, x >= 0
%==============================================
% where
H = [4 0 0; 0 1 -1; 0 -1 1];
p = [-8 -6 -6]';
A = [1 1 1];
b = 3;
% using the initial point initial point (x0,lmd0,mu0) where
x0 = [1 2 2]';
lmd0 = -1;
mu0 = [0.2 0.2 0.2]';
rho = 3;
epsi = 1e-6;
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