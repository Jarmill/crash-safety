function [box] = poly_bounding_box(A, b)
%POLY_BOUNDING_BOX Find a bounding box of the convex polytope Aw <= b

%by linear programming of polytope

[m, n] = size(A);

box = zeros(n, 2);

%linear programming options
options = optimset;
options = optimset(options,'Display', 'off');

for i = 1:n
    %objectives for dimension i
    f_pos = full(sparse(i, 1, -1, n, 1));
    f_neg = full(sparse(i, 1, 1, n, 1));
    
    
    %linprog is a minimization objective
    %be careful of signs
    c_pos = linprog(f_pos,A,b,[],[],[],[],[],options);
    c_neg= linprog(f_neg,A,b,[],[],[],[],[],options);
    
    bound_pos = -f_pos'*c_pos;
    bound_neg = f_neg'*c_neg;
    
    %record the box
    box(i, :) = [bound_neg, bound_pos];
    
end


end

