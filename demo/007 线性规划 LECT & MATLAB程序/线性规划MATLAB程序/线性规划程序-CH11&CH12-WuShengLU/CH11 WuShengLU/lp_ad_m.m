% Program: lp_ad_m.m
% Title: Simplex algorithm for alternative-form LP
% problems with degenerate vertices, that does not require
% an initial point as an input.
% Description: implements Algorithm 11.2 for alternative-form 
% LP problems with degenerate vertices, that does not require
% an initial point as an input.
% Theory: See Practical Optimization Secs. 11.3.1.2 and 11.2.3.4.
% Input:  
%    A, b, c - input data matrices
% Output:        
%    xs - a vertex minimizer.
%    fs - objective function at xs
%    k - number of iteration at convergence
% Example:
% Find a minimizer of the alternative-form LP problem
%    minimize c'*x 
%    subject to A*x >= b
% where
%    A = [-1 -1 -1; -0.5 -2 -1; eye(3)]
%    b = [-1 -1 0 0 0]'
%    c = [0 -1 0]'
% Solution:
% Execute the command
%    [xs,fs,k] = lp_ad_m(A,b,c)
% ====================================================
function [xs,fs,k] = lp_ad_m(A,b,c)
disp('    ')
disp('Program lp_ad_m.m')
k = 0;
[p,n] = size(A);
ep = ones(p,1);
zn = zeros(n,1);
Aw = [A ep];
Ae = [Aw; [zn' 1]];
phi0 = max([0; b]);
xe0 = [zn; phi0];
be = [b; 0];
ce = [zn; 1];
[xes,fes,ke] = lp_ad(Ae,be,ce,xe0);
if fes == 0,
    x0 = xes(1:n);
else
    disp('Feasible region is empty.')
    xs = [];
    fs = [];
return
end
x = find_v(A,b,x0);
r = A*x-b;
I0 = (1:1:p)';
if r >= -1e-10,
else
disp('Point x0 is not feasible.')
xs = [];
fs = [];
return
end
r(find(r < 1e-10)) = 0;
act_set = find(r == 0);
Lact = length(act_set);
if Lact >= n,
   [Qs,Rs,E] = qr(A(act_set,:)');
   I1 = sum(E(:,1:n)')';
   I2 = I1.*act_set;
   I3 = find(I2 ~= 0);
   act_set = I2(I3);
   Aa = A(act_set,:);
else
   disp('Point x is not a vertex.')
   xs = [];
   fs = [];
   return
end
if norm(x - x0) > 1e-12,
   k = 1;
end
while 1,   
  mu = inv(Aa')*c;
   if mu > -1e-10,
      break
   end
  I4 = find(mu <= -1e-10);
  l = I4(1);
  iAa = inv(Aa);
  d = iAa(:,l);
  Ad = A*d;
  I5 = setdiff(I0,act_set);
  I6 = find(Ad(I5) <= -1e-10);
  I = I5(I6);
   if isempty(I),
    disp('The objective function unbounded from below.')
    xs = [];
    fs = [];
    return
   end
  [alpha, ii] = min(r(I)./(-A(I,:)*d));
  act_set(l) = I(ii);
  act_set = sort(act_set);
  Aa = A(act_set,:);
  x = x + alpha*d;
  r = A*x - b;
  k = k + 1;
end
xs = x;
fs = c'*xs;