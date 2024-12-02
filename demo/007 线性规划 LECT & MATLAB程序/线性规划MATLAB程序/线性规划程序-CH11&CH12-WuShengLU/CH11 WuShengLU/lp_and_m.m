% Program: lp_and_m.m
% Title: Simplex algorithm for alternative-form LP
% problems with non-degenerate vertices, that does
% not require an initial point as an input.
% Description: implements Algorithm 11.1 for alternative-form 
% LP problems with non-degenerate vertices, that does not 
% require an initial point as an input.
% Theory: See Practical Optimization Secs. 11.3.1.1 and 11.2.3.4.
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
%    A = [1 0; -1 0; 0 1; -1 -1; -1 -2]
%    b = [0 -2 0 -3.5 -6]'
%    c = [-1 -4]'
% Solution:
% Execute the command
%    [xs,fs,k] = lp_and_m(A,b,c)
% ======================================================
function [xs,fs,k] = lp_and_m(A,b,c)
disp('    ')
disp('Program lp_and_m.m')
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
[xes,fes,ke] = lp_and(Ae,be,ce,xe0);
if fes == 0,
    x0 = xes(1:n);
else
    disp('The feasible region is empty.')
    xs = [];
    fs = [];
return
end
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