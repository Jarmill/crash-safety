function X_cell = prep_space_cell(X)
    %PREP_SPACE_CELL: Split up the support set X into cells
    if iscell(X)
        %X is already a cell, nothing to be done
        X_cell = cell(length(X), 1);
        for i = 1:length(X)
            X_cell{i} = fill_constraint(X{i});
        end
    else
        if isstruct(X)
            %wrap up the structure in a cell
            X_cell = {X};
        else
            %X is supported at discrete points
            %Each column of X is a possible origin datapoint
            if size(X, 2) == 1
                X_cell = {X};
            else
                %each cell element is a column of the original array 
                Xt = X';
                Xt_cell = mat2cell(Xt, ones(size(X, 2), 1));
                X_cell = cellfun(@(x) x', Xt_cell, 'UniformOutput', false);
            end
        end
    end
end