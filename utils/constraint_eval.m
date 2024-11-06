function [cons_eval, handle] = constraint_eval(X, vars, pt, eq_tol)
%CONSTRAINT_EVAL Evaluate the constraints X.ineq>=0, X.eq==0 at point pt

if nargin < 4
    eq_tol = 1e-8;
end

% ineq = replace(X.ineq, vars, pt);
% eq   = replace(X.eq, vars, pt);

ineq_func = polyval_func(X.ineq, vars);
ineq = ineq_func(pt);

eq_func = polyval_func(X.eq, vars);
eq = eq_func(pt);


ineq_eval = ineq >= 0;
eq_eval = abs(eq) <= eq_tol;

cons_eval = all([ineq_eval; eq_eval]);

handle = @(x) all([ineq_func(x) >= 0; abs(eq_func(x)) <= eq_tol]);

end

