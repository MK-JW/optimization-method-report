% Program: kelley_ie.m
% Title: Kelley's cutting-plane algorithm for the CP 
% problem in Eq, (13.70).
% Description: Implements Algorithm 13.6 for the CP 
% problem in Eq. (13.70).
% Theory: See Practical Optimization Sec. 13.5.
% Input:     
%         cname -- name of MATLAB function that evaluates
%                  functions cj(x)for j = 1, 2, ..., q in 
%                  Eq.(13.70b).
%         hname -- name of MATLAB function that evaluates
%                  gradients of -cj(x) for j = 1,2,...,q.
%            X0 -- matrix whose columns are a set of 
%                  feasible points
%          epsi -- convergence tolerance
% Output:   
%            xs -- solution vector
%            fs -- value of objective function at xs
%             k -- number of iterations at convergence
% Example:
% Apply Algorithm 13.6 to solve the QP problem
% minimize 0.5*x'*H*x +x'*p
% subject to 
% -x1 + x2 + 2 >= 0
% -x1 - x2 + 6 >= 0
% x1 >= 0
% x2 >= 0
% where
% H = [4 1; 1 1]
% p = [-1 2]'
% Solution:
% First, convert the problem at hand to the form of
% Eq. (13.70) where
% x = [x1 x2 L]'
% c1(x) = L - (0.5*x'*H*x +x'*p)
% c2(x) = -x1 + x2 + 2
% c3(x) = -x1 - x2 + 6
% c4(x) = x1
% c5(x) = x2
% Function  kelley_ie requires two MATLAB functions cex1.m 
% and hex1.m as follows:
% ====================
% function c = cex1(x)
% x1 = x(1);
% x2 = x(2);
% L = x(3);
% xx = x(1:2);
% H = [4 1; 1 1];
% p = [-1 2]';
% c1 = L - (0.5*xx'*H*xx +xx'*p);
% c2 = -x1 + x2 + 2;
% c3 = -x1 - x2 + 6;
% c4 = x1;
% c5 = x2;
% c = [c1 c2 c3 c4 c5]';
% ======================
% function h = hex1(x)
% xx = x(1:2);
% H = [4 1; 1 1];
% p = [-1 2]';
% h1 = [H*xx+p; -1];
% h2 = [1 -1 0]';
% h3 = [1 1 0]';
% h4 = [-1 0 0]';
% h5 = [0 -1 0]';
% h = [h1 h2 h3 h4 h5];
% =======================
% These m-files can be found in the same folder as m-file kelley_ie.m.
% The solution can be obtained by executing the following commands:
% X0 = [1 1 2 2; 1 2 1 2; 5 10 11 18]
% epsi = 1e-7
% [xs,fs,k]= kelley_ie('cex1','hex1',X0,epsi)
% =================================================
function [xs,fs,k] = kelley_ie(cname,hname,X0,epsi)
disp(' ')
disp('Program kelley_ie.m')
[n,K] = size(X0);
Ak = [];
bk = [];
c = [zeros(n-1,1);1];
x0 = X0(:,1);
for i = 1:K,
    xi = X0(:,i);
    ci = feval(cname,xi);
    hi = feval(hname,xi);
    Ai = -hi';
    Ak = [Ak; Ai];
    bk = [bk; Ai*xi-ci];
end
k = 1;
[xw,fsw,kw] = lp_ad(Ak,bk,c,x0);
cw = feval(cname,xw);
[yw,iw] = min(cw);
while yw < -epsi,
    cjs = cw(iw);
    h1 = feval(hname,xw);
    hjs = h1(:,iw);
    Ak = [Ak; -hjs'];
    bk = [bk; -hjs'*xw-cjs];
    [xw,fsw,kw] = lp_ad(Ak,bk,c,x0);
    cw = feval(cname,xw);
    [yw,iw] = min(cw);
    k = k + 1;
end
xs = xw;
fs = xs(end); 