% Piecewise linear minimization via Ellipsoid algorithm with constraints 
% EE364b example
% 
% PWL minimization problem: min max_i=1..m (a_i'x + b_i)
%    subject to  -0.1 <= x <= 0.1
%
% number of variables n and linear functions m
n = 20; m = 100;

%********************************************************************
% generate an example
%********************************************************************
randn('state',1)
A = randn(m,n);
b = randn(m,1);
C = [eye(n,n);-eye(n,n)];   % constraints  -0.1 <= xi <= 0.1
c = 0.1*ones(2*n,1);        % constraints  C*x <= c
%********************************************************************
% compute pwl optimal point using CVX
%********************************************************************
cvx_begin
    variable x(n)  
    minimize max(A*x + b)
    subject to
       x <= 0.1
       x >= -0.1
cvx_end
f_star = cvx_optval;

%********************************************************************
% ellipsoid method with constraints
%********************************************************************
eps = 0.01;
% number of iterations 
niter = 2000;
% initial ellipsoid
P = eye(n);
% initial point 
x = zeros(n,1); 
U = [+inf]; L = [-inf];
Lm = [-inf];

f_save = []; f_best = []; 



for iter = 1:niter 
    % find active functions at current x
    [f, idx] = max(A*x + b);        
    % subgradient at current x
    g = A(idx(1),:)';
	% convergence test
	if sqrt(g'*P*g) <= eps && max(C*x-c) <= eps
		disp('converged')
		break
	end
    % save current and best function values 
    f_save = [f_save f];
    f_best = [f_best min(f_save)];

	% update upper / lower bounds.
	U(iter+1) = min(U(iter), max(A*x + b));
	L(iter+1) = max(L(iter), max(A*x + b) - sqrt(g'*P*g));
	Lm(iter+1) = max(A*x + b) - sqrt(g'*P*g);

	h = f - f_best(end);

    [fc, idc] = max(C*x-c); % Find the maximum value of the constraint function
    if fc < 0 
        % x is feasible
        % update with deep cut.
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
    else 
        % fi(x) is infeasible for some i
        % the most violated constraint is f_idc
        g = C(idc(1),:)';
        % update with deep cut.
        h=C(idc(1),:)*x-c(idc(1));
        alpha = h/sqrt(g'*P*g);
		gt = g/sqrt(g'*P*g);
		x = x - (1+alpha*n)/(n+1)*P*gt;
		P = n^2*(1-alpha^2)/(n^2-1)*(P-2*(1+alpha*n)/((n+1)*(1+alpha))*P*(gt*gt')*P);
    end
end

% trim U, L. Unfortunate hack.
U = U(2:end);
L = L(2:end);
Lm = Lm(2:end);

x
