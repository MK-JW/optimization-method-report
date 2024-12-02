function [x_optimal, f_optimal, exit_flag, output] = Conjugate_gradient_optimization()
    % Define the initial point and tolerance
    x_initial = [4; -4; 4; -4];
    tolerance = 1e-6;

    % Use fminunc for optimization (Conjugate Gradient method)
   % options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
   % options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
   options = optimset('GradObj', 'on', 'MaxIter', 1000, 'TolFun', tolerance, 'TolX', tolerance);
    problem.options = options;
    problem.x0 = [4; -4; 4; -4];
    problem.objective = @f_testwithgrad;
    problem.solver = 'fminunc';
    [x_optimal, f_optimal, exit_flag, output]  = fminunc(problem)   
    
    % Display results
    disp(['Optimal x: ', mat2str(x_optimal)]);
    disp(['Optimal f: ', num2str(f_optimal)]);
    disp(['Exit flag: ', num2str(exit_flag)]);
    disp(['Number of iterations: ', num2str(output.iterations)]);
end

%function f = f_test3(x)
  %  x1 = x(1);
  %  x2 = x(2);
  %  x3 = x(3);
  %  x4 = x(4);
 % syms x1 x2 x3 x4
 %   f = 100 * (x1^2 - x2)^2 + (1 - x1)^2 + (1 - x3^2)^2 + 10 * ((x2 - 1)^2 + (x4 - 1)^2) + 18 * (x2 - 1) * (x4 - 1);
 %   g=jacobian(f)
%    f = 100 * (x(1)^2 - x(2))^2 + (1 - x(1))^2 + (1 - x(3)^2)^2 + 10 * ((x(2) - 1)^2 + (x(4) - 1)^2) + 18 * (x(2) - 1) * (x(4) - 1);
%end

%[x_optimal, f_optimal, exit_flag, output] = Conjugate_gradient_optimization()

function [f,g] = f_testwithgrad(x)
% Calculate objective f
f = 100 * (x(1)^2 - x(2))^2 + (1 - x(1))^2 + (1 - x(3)^2)^2 + 10 * ((x(2) - 1)^2 + (x(4) - 1)^2) + 18 * (x(2) - 1) * (x(4) - 1);

if nargout > 1 % gradient required
    g = [ 2*x(1) - 400*x(1)*(- x(1)^2 + x(2)) - 2;
        - 200*x(1)^2 + 220*x(2) + 18*x(4) - 38;
        4*x(3)*(x(3)^2 - 1);
        18*x(2) + 20*x(4) - 38];
end
end
