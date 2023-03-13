%what is the relationship between the number of nonredundant faces on a
%polytope and dilation of the rhs constraints?

%increasing z could actually increase the number of faces in the polytope.
%this is not the result that I was expecting.
%

load('crash_poly_flow.mat');

[m, L] = size(W.A);

w = sdpvar(L, 1);
z = sdpvar(1, 1);

%find the first nontrivial parameter value

cons = [z>=0; W.A*w <= (z + W.b)];

opts = sdpsettings('solver', 'mosek', 'verbose', 0);
sol = optimize(cons, z, opts);

z_min = value(z)

[A_min, b_min, non_min, num_min] = nontrivial_constraints(W.A, W.b+z_min);
[A_half, b_half, non_half, num_half] = nontrivial_constraints(W.A, W.b+0.5);
[A_opt, b_opt, non_opt, num_opt] = nontrivial_constraints(W.A, W.b+0.55);
[A_all, b_all, non_all, num_all] = nontrivial_constraints(W.A, W.b+1);
[A_double, b_double, non_double, num_double] = nontrivial_constraints(W.A, W.b+2);