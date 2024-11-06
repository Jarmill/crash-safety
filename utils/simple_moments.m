function [mom_out] = simple_moments(var_count, order)
%SIMPLE_MOMENTS monomial powers and moments without symmetry
%Inputs:
%  var_count: the number of variables
%
%
% order:    degree of relaxation, monomial powers
%           moment matrix will have moments of up to degree d=2*order
%
%Outputs:   struct mom_out with fields:
%   monom_int:      monom_int: monomials that are on the interior of the
%                   moment matrices
%   index:          Don't know what this does
%   M:              indices for moments (vectorized, not a matrix)
%   N:              Number of moments on the interior of the moment matrix
%   order:          order of relaxation
%% preliminary counting and indexing

d = 2*order;

mom_out = struct;

%% Form monomial basis 

%find the monomial basis
[basis, ~] = momentPowers(0, var_count, order);

%% Find monomial powers inside the matrix

m = size(basis, 1);
temp = bsxfun(@plus,kron(ones(m, 1), basis), kron(basis, ones(m, 1)));

%monomials indexed by moment matrix
monom_int = unique(temp, 'rows');
index = 1:size(monom_int, 1);

% index_order2 = any(monom_int(:, 1:sum(var_count(1:2))) == 2, 2);

%vectorized moment matrix
M = get_vec_ind_func(monom_int,temp,index,0)';

% ZZ_even = unique(temp_even, 'rows');
mom_out.monom_int = monom_int;
mom_out.index = index;
% mom_out.index_order2 = index_order2; 
% mom_out.M_even = M_even;
mom_out.M = M;

%number of monomials
mom_out.N_int = index(end);
mom_out.N= m;
% mom_out.N_odd = m_odd;
mom_out.order = order;

end

