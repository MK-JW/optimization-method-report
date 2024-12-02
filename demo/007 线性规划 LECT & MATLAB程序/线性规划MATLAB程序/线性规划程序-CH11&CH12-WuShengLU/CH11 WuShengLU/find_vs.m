% Program: find_vs.m
% Description: Finds a vertex of the feasible region
% defined by A*x = b, and x >= 0.
% Theory: See Practical Optimization Secs.11.2.3.3
%    and 11.2.3.4.
% Input:  
%   A, b - input data matrices
%   x0 - a feasible initial point
% Output:   
%   x - a vertex of the feasible region
% Example:
% Find a vertex of the feasible region described
% by A*x = b and x >= 0 where 
%    A = [1 1 1]
%    b = 1
% The initial point is given by
%    x0 = [1/3 1/3 1/3]'
% Solution:
% Execute the command
%    x = find_vs(A,b,x0)
% ==================================================
function x = find_vs(A,b,x0)
disp('    ')
disp('Program find_vs.m')
[p, n] = size(A);
[Q,R] = qr(A');
R = R(1:p,:);
Q1 = Q(:,1:p);
Q2 = Q(:,p+1:n);
xs = Q1*inv(R')*b;
phi0 = Q2'*(x0-xs);
phi = find_v(Q2,-xs,phi0);
x = Q2*phi + xs;