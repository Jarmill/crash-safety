function [supcon_out] = constraint_convert(Cons, var, varnew)
%CONSTRAINT_CONVERT Convert constraints from yalmip form (struct with ineq
%and eq) to gloptipoly @supcon

f_ineq = polyval_func(Cons.ineq, var);
f_eq = polyval_func(Cons.eq, var);

supcon_out = [f_ineq(varnew)>=0; f_eq(varnew) == 0];


end

