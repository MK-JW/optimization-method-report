% Program: qp_e.m
% Title: QR-decomposition based algorithm for
% convex QP problems with equality constraints.
% Description: Implements a QR-decomposition
% algorithm for convex QP problems with 
% inequality constraints.
% Theory: See Practical Optimization Sec. 13.2.
% Input:  
%     H -- positive semidefinite Hessian matrix
%     A -- full row-rank constraint matrix A
%  p, b -- input vectors
% Output:        
%          xs -- solution vector
% Example:
% Find a minimizer of the QP problem
% minimize 0.5*x'*H*x +x'*p
% subject to A*x = b
% where
% H = [1 0 0; 0 1 0; 0 0 0.01]
% p = [2 1 -1]'
% A = [0 1 1]
% b = 1
% Solution:
% Execute the following commands:
% H = [1 0 0; 0 1 0; 0 0 0.01]
% p = [2 1 -1]'
% A = [0 1 1]
% b = 1
% xs = qp_e(H,p,A,b)
% =============================
function xs = qp_e(H,p,A,b)
% disp(' ')
% disp('Program qp_e.m')
if min(size(A))== 0,
 xs = -pinv(H)*p;
else
 n = length(p);
 p1 = size(A)*[1 0]';
 [Q,R] = qr(A');
 R = R(1:p1,:);
 Q1 = Q(:,1:p1);
 Q2 = Q(:,p1+1:n);
 xs = Q1*inv(R')*b;
 Hh = Q2'*H*Q2;
 ph = Q2'*(H*xs+p);
 phi_s = -inv(Hh)*ph;
 xs = Q2*phi_s+xs;
end