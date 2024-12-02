% Program: lp_and.m
% Title: Simplex algorithm for alternative-form LP
% problems with non-degenerate vertices.
% Description: Implements Algorithm 11.1 for alternative-form 
% LP problems with non-degenerate vertices.
% Theory: See Practical Optimization Sec. 11.3.1.1.
% Input:  
%    A, b, c - input data matrices
%    x0 - a feasible initial point
% Output:        
%    xs - a vertex minimizer.
%    fs - objective function at xs
%    k - number of iteration at convergence
% Example:
% Find a minimizer of the alternative-form LP problem
%    minimize c'*x 
%    subject to A*x >= b
% where
%    A = [1 0; -1 0; 0 1; -1 -1; -1 -2]
%    b = [0 -2 0 -3.5 -6]'
%    c = [-1 -4]'
% The initial point is given by
%    x0 = [1 1]'
% Solution:
% Execute the command
%    [xs,fs,k] = lp_and(A,b,c,x0)
% =========================================================
function [xs,fs,k] = lp_and(A,b,c,x0)
disp('    ')
disp('Program lp_and.m')
k = 0;
[p, n] = size(A);
I0 = (1:1:p)';
x = find_v(A,b,x0);
r = A*x-b;
if r >= -1e-10,
else
  disp('Point x is not feasible.')
  xs = [];
  fs = [];
return
end
e = eye(n);
r(find(r < 1e-10)) = 0;
act_set = find(r == 0);
Aa = A(act_set,:); 
if rank(Aa) == n,
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
 mu = Aa'\c;
 if mu > -1e-10,
    break
 end
 I1 = find(mu < -1e-10);
 l = I1(1);
 d = Aa\e(:,l);
 Ad = A*d;
 I2 = setdiff(I0,act_set);
 I3 = find(Ad(I2) <= -1e-10);
 I = I2(I3);
  if isempty(I),
   disp('The objective function is unbounded from below;')
   disp('hence no solution exists.')
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