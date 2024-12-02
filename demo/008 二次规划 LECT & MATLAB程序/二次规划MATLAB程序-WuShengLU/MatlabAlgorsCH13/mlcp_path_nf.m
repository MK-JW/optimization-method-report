% Program: mlcp_path_nf.m
% Title: Nonfeasible-initialization interior-point 
% algorithm for mixed LCP problems
% Description: Implements Algorithm 13.4 for mixed 
% LCP problems
% Theory: See Practical Optimization Sec. 13.4.4.
% Input:   
% {K11,K12,K21,K22) -- data matrices to form matrix K
%           {q1,q2} -- data vectors to form vector q
%     (x0,lmd0,mu0) -- initial point with x0 > 0 and mu0 > 0
%               rho -- parameter involved in determining 
%                      the value of tau, rhou is required 
%                      to be no less than the square root of n
%                ep -- tolerance for duality gap
% Output:   
%        (x,mu,lmd) -- solution vectors
%                 k -- number of iterations at convergence
% Example:
% Apply Algorithm 13.4 to solve the QP problem
% minimize 0.5*x'*H*x +x'*p
% subject to A*x = b, x >= 0
% where
% H = [4 0 0; 0 1 -1; 0 -1 1]
% p = [-8 -6 -6]'
% A = [1 1 1]
% b = 3
% Solution:
% Execute the following commands:
% H = [4 0 0; 0 1 -1; 0 -1 1]
% p = [-8 -6 -6]'
% A = [1 1 1]
% b = 3
% K11 = H;
% K12 = -A';
% K21 = A;
% K22 = 0;
% q1 = p;
% q2 = -b;
% x0 = [1 2 2]'
% lmd0 = -1
% mu0 = [0.2 0.2 0.2]'
% rho = 3
% epsi = 1e-6
% [x,mu,lmd,k]= mlcp_path_nf(K11,K12,K21,K22,q1,q2,x0,mu0,lmd0,rho,epsi)
% ===================================================
function [x,mu,lmd,k]= mlcp_path_nf(K11,K12,K21,K22,q1,q2,x0,...
    mu0,lmd0,rho,epsi)
disp(' ')
disp('Program mlcp_path_nf.m')
% Data initialization.
n = length(x0);
p = length(lmd0);
rho_n = n + rho;
a_max = 1-1e-6;
x = x0(:);
lmd = lmd0(:);
mu = mu0(:);
xm = x.*mu;
gap = sum(xm);
I = eye(n);
Z = zeros(p,n);
e = ones(n,1);
k = 0;
% Iteration begins.
while gap > epsi,
  % Compute tau, r1, and r2.
  tau = gap/rho_n;
  r1 = K11*x + K12*lmd - mu + q1;
  r2 = -K21*x - K22*lmd -q2;
  M = diag(mu);
  X = diag(x);
  % Calculate d_x,d_lmd, and d_mu.
  C = [-K11 -K12 I; K21 K22 Z; M Z' X];
  rr = [r1; r2; tau*e-xm];
  dw = inv(C)*rr;
  d_x = dw(1:n);
  d_lmd = dw((n+1):(n+p));
  d_mu = dw((n+p+1):(2*n+p));
  % Calculate step size.
  ind = find(d_x < 0);
  a_p = min(x(ind)./(-d_x(ind)));
  ind = find(d_mu < 0);
  a_d = min(mu(ind)./(-d_mu(ind)));
  a_k = a_max*min([a_p a_d]);
  % Form new iterate.
  x = x + a_k*d_x;
  lmd = lmd + a_k*d_lmd;
  mu = mu + a_k*d_mu;
  % Compute duality gap and check convergence.
  xm = x.*mu;
  gap = sum(xm);
  % Update iteration index, etc.
  k = k + 1;
end