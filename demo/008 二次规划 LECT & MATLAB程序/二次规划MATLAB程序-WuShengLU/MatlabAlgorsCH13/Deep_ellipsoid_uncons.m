function [xs,fs,f_save,f_best,Lm,L,U] = Deep_ellipsoid_uncons(fname,gname,x0,A0,DEEPCUT,eps)
disp(' ')
disp('Program Deep_ellipsoid_uncons.m')
format long
%
%********************************************************************
% ellipsoid method
%********************************************************************
%eps = 0.01;
% number of iterations 
niter = 2000;
% initial ellipsoid
if A0 == []
    P = eye(n);
else
    P=A0;
end
% initial point 
x = x0; 
U = [+inf]; L = [-inf];
Lm = [-inf];

f_save = []; f_best = []; 
for iter = 1:niter 
    % find function value at current x
    f = feval(fname,x);
    % subgradient at current x
    g = feval(gname,x);
	% convergence test
	if sqrt(g'*P*g) <= eps
        xs=x;fs=f;
		disp('converged')
		break
	end
    % save current and best function values 
    f_save = [f_save f];
    f_best = [f_best min(f_save)];
   	% update upper / lower bounds.
    %f = feval(fname,x);
	U(iter+1) = min(U(iter), f);
	L(iter+1) = max(L(iter), f - sqrt(g'*P*g));
	Lm(iter+1) = f - sqrt(g'*P*g);

	h = f - f_best(end);

	if DEEPCUT && h > 0
		% update with deep cut.
		alpha = h/sqrt(g'*P*g);
		gt = g/sqrt(g'*P*g);
		x = x - (1+alpha*n)/(n+1)*P*gt;
		P = n^2*(1-alpha^2)/(n^2-1)*(P-2*(1+alpha*n)/((n+1)*(1+alpha))*P*(gt*gt')*P);
	else
		% update ellipsoid with shallow cut.
		gt = g/sqrt(g'*P*g);
		x = x - 1/(n+1)*P*gt;
		P = n^2/(n^2-1)*(P-2/(n+1)*P*(gt*gt')*P);
    end
    
end

% trim U, L. Unfortunate hack.
U = U(2:end);
L = L(2:end);
Lm = Lm(2:end);
end
