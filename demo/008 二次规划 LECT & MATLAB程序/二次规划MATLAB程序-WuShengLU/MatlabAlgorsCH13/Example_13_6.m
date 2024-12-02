function Example_13_6
close;clear;clc;
% =====Initial ==================================================
epsi=1e-7;

%x0=[4,-1]'
%A0=[36,0;0,9]
%[xs,fs,k] = ellipsoid_uncons_CP('fex','gex',x0,A0,epsi)

x0=[31,11]'
A0=[1600,0;0,256]
[xs,fs,k] = ellipsoid_uncons_CP('fex','gex',x0,A0,epsi)

%x0=[61,-29]'
%A0=[6400,0;0,2304]
%[xs,fs,k] = ellipsoid_uncons_CP('fex','gex',x0,A0,epsi)

end

% Program: ellipsoid_ie.m
% Title: Ellipsoid method for unconstrained CP problem
% Description: Implements Algorithm 13.7 for CP problem
% Theory: See Practical Optimization Sec. 13.6.
% Input: 
%         fname -- name of MATLAB function that evaluates
%                  objective function f(x)
%         gname -- name of MATLAB function that evaluates
%                  the gradient of objective function f(x)
%      (x0, A0) -- input data that characterize an initial 
%                  ellipsoid {x: (x-x0)'*inv(A0)*(x-x0)} 
%                  that contains a minimizer
%          epsi -- convergence tolerance
% Output:   
%            xs -- solution vector
%            fs -- value of objective function at xs
%             k -- number of iterations at convergence
% Example:13.6
% Apply Algorithm 13.7 to solve the unconstrainted QP problem
% minimize f(x)
    function [xs,fs,k] = ellipsoid_uncons_CP(fname,gname,x0,A0,epsi)
disp(' ')
disp('Program ellipsoid_ie.m')
format long
n = length(x0);
x = x0;
A = A0;
gak = 1;
k = 0;

xk_best=x0;
fk_best =feval(fname,xk_best);

gk = feval(gname,x);
gak = sqrt(gk'*A*gk);
%--------------figure------------------------
figure
plot(0,feval(fname,x0),'or')
%--------------------------------------------
while gak > epsi
      gkt = gk/gak;
      x = x - A*gkt/(n+1);
      Az = A*gkt;
      A = n^2*(A - 2*(Az*Az')/(n+1))/(n^2-1);
      fk = feval(fname,x);
      if fk < fk_best
          fk_best=fk;
          xk_best=x;
      end
      k = k + 1; 
      gk = feval(gname,x);     
      gak = sqrt(gk'*A*gk);
   
%--------------figure------------------------
hold on
plot(k,fk_best,'ob')
%--------------------------------------------
 end
%xs = x;
%fs = feval(fname,xs);
xs=xk_best;
fs=fk_best;
end

% ============= Evaluate a subgradient at x =============
% f(x) = (x1 − 5*x2 + 4)^2 + (7*x1 + 11*x2 − 18)^4
 function g = gex(x)
 x1 = x(1)-5*x(2)+4; 
 x2= (7*x(1) + 11*x(2)-18)^3;
 g1=  2*x1 + 28*x2;
 g2=-10*x2 + 44*x2;
 g=[g1;g2];
 end
% ========================================================

% ========================================================
% f(x) = (x1−5*x2+4)^2+(7*x1+11*x2−18)^4
function f = fex(x)
f = (x(1)-5*x(2)+4)^2 + (7*x(1)+11*x(2)-18)^4;
end
% =====================
% 