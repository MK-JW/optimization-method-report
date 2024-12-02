% Program: ellipsoid_ie.m
% Title: Ellipsoid method for constrained CP problem
% in Eq, (13.79).
% Description: Implements Algorithm 13.8 for CP 
% problem in Eq. (13.79).
% Theory: See Practical Optimization Sec. 13.6.
% Input: 
%         fname -- name of MATLAB function that evaluates
%                  objective function f(x)
%         gname -- name of MATLAB function that evaluates
%                  the gradient of objective function f(x)
%         cname -- name of MATLAB function that evaluates
%                  constraint functions cj(x)for j = 1, 
%                  2, ..., q in Eq.(13.79b).
%         hname -- name of MATLAB function that evaluates
%                  the gradients of -cj(x) for j = 1,2,...,q.
%      (x0, A0) -- input data that characterize an initial 
%                  ellipsoid {x: (x-x0)'*inv(A0)*(x-x0)} 
%                  that contains a minimizer
%          epsi -- convergence tolerance
% Output:   
%            xs -- solution vector
%            fs -- value of objective function at xs
%             k -- number of iterations at convergence
% Example:
% Apply Algorithm 13.8 to solve the QP problem
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
% Function ellipsoid_ie, requires four MATLAB function m-files  
%  fex.m, gex.m, cex.m and hex.m as follows:
% =====================
% function f = fex(x)
% H = [4 1; 1 1];
% p = [-1 2]';
% f = 0.5*x'*H*x + x'*p;
% =====================
% function g = gex(x)
% H = [4 1; 1 1];
% p = [-1 2]';
% g = H*x + p;
% ====================
% function c = cex(x)
% x1 = x(1);
% x2 = x(2);
% c1 = -x1 + x2 + 2;
% c2 = -x1 - x2 + 6;
% c3 = x1;
% c4 = x2;
% c = [c1 c2 c3 c4]';
% ======================
% function h = hex(x)
% h1 = [1 -1]';
% h2 = [1 1]';
% h3 = [-1 0]';
% h4 = [0 -1]';
% h = [h1 h2 h3 h4];
% =======================
% These m-files are included in the same folder as m-file
% ellipsoid_ie.m. 
% Next, let
% x0 = [1.5 3]'
% A0 = 8*eye(2,2)
% epsi = 1e-7
% and execute the command
% [xs,fs,k]= ellipsoid_ie('fex','gex','cex','hex',x0,A0,epsi)
% =================================================
function [xs,fs,k] = ellipsoid_ie(fname,gname,cname,hname,...
x0,A0,epsi)
disp(' ')
disp('Program ellipsoid_ie.m')
n = length(x0);
x = x0;
A = A0;
gak = 1;
k = 0;
while gak > epsi,
  c = feval(cname,x);
  [y,ind] = min(c);
  if y >= 0,
     gk = feval(gname,x);
     gak = sqrt(gk'*A*gk);
       if gak < epsi,
          xs = x;
       break
       else
          gkt = gk/gak;
          x = x - A*gkt/(n+1);
          Az = A*gkt;
          A = n^2*(A - 2*Az*Az'/(n+1))/(n^2-1);
          k = k + 1;
       end
  else
      h = feval(hname,x);
      hk = h(:,ind);
      hak  = sqrt(hk'*A*hk);
      gkt = hk/hak;
      zz = y + hak;
      if zz < 0,
         disp('No solution was found.')
         break
      else
         x = x - A*gkt/(n+1);
         Az = A*gkt;
         A = n^2*(A - 2*Az*Az'/(n+1))/(n^2-1);
         k = k + 1;
      end
   end
end
xs = x;
fs = feval(fname,xs);