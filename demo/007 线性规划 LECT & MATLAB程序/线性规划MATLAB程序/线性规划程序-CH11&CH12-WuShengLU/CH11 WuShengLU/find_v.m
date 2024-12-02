% Program: find_v.m
% Description: Finds a vertex of the feasible region
% defined by A*x >= b.
% Theory: See Practical Optimization Secs.11.2.3.3 
% and 11.2.3.4.
% Input:  
%   A, b - input data matrices
%   x0 - a feasible initial point
% Output:    
%   x - a vertex of the feasible region
% Example:
% Find a vertex of the feasible region described
% by A*x >= b where 
%    A = [1 0; -1 0; 0 1; -1 -1; -1 -2]
%    b = [0 -2 0 -3.5 -6]'
% The initial point is given by
%    x0 = [1 1]'
% Solution:
% Execute the command
% x = find_v(A,b,x0)
% ====================================================
function x = find_v(A,b,x0)
disp('    ')
disp('Program find_v.m')
[p,n] = size(A);
x = x0;
 if A*x0-b >= -1e-14,
 else
  disp('Point x0 is not feasible.')
  x = [];
  return
 end
 r = A*x-b;
 r(find(r < 1e-14)) = 0;
 act_set = find(r == 0);
 Lact = length(act_set);
 if Lact ~= 0,
   Aa = A(act_set,:);
 else
   Aa = [];
 end
while 1,
   ra = rank(Aa);
   dr = Lact-ra;
   if dr > 0,
      [Q1,R1,E] = qr(Aa');
      I1 = sum(E(:,1:ra)')';
      I2 = I1.*act_set;
      I3 = find(I2 ~= 0);
      act_set = I2(I3);
      Aa = A(act_set,:);
   end
   if ra == n,
      break
   end
   inact_set = find(r > 0);
   Lin = length(inact_set);
   wi = 1;
   si = 0;
   for i = 1:Lin,
    ap = [Aa; A(inact_set(i),:)];
    [Q2,R2] = qr(ap');
     cp = R2(ra+1,ra+1);
      if abs(cp) > si,
       wi = i;
       si = abs(cp);
      end
    end
   ap = [Aa; A(inact_set(wi),:)];
   d = pinv(ap)*[zeros(size(Aa,1),1); -1];
   I = find((A*x-b > 0) & (A*d < 0));
    if length(I) ~= 0,
     [alpha, ix] = min((A(I,:)*x-b(I))./(-A(I,:)*d));
     x = x +alpha*d;
     act_set = [act_set;I(ix)];
     Lact = length(act_set);
     Aa = A(act_set,:);
     r = A*x-b;
    end
end