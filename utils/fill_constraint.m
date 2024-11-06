function [con_out] = fill_constraint(con_in)
%EMPTY_CONSTRAINT
%con_in should be a struct with fields 'ineq' and 'eq'
%add them if they are missing
    con_out = con_in;
    if isempty(con_in) 
        con_out = struct('ineq', [], 'eq', []);        
    elseif ~isa(con_in, 'struct')
        con_out = struct('ineq', con_in, 'eq', []);
    else    
        if ~isfield(con_in, 'ineq')
            con_out.ineq = [];
        end

        if ~isfield(con_in, 'eq')
            con_out.eq = [];
        end
    end
end
