function [box_out, box_center, box_half] = box_process(nx, box)
%BOX_PROCESS Process the box constraint in box:
%nx: number of states in system (length(x))
%box: input box constraint
    %box is scalar B:           -B     <= x_i <= B
    %box is [B_i]:              -Bi    <= x_i <= B_i
    %box is [Bmin, Bmax]:       Bmin   <= x_i <= Bmax    
    %box is [Bmin_i, Bmax_i]    Bmin_i <= x_i <= Bmax_i
%box_process turns box into the 4th form [Bmin_i, Bmax_i]
%
%Output:
%   box_out:    box output, as mentioned above
%   box_center: center of box in each coordinate
%   box_half:   half-length of box from center


%box_out = zeros(nx, 2);

assert(size(box, 1)==nx || size(box, 1)==1);

if size(box, 1) == 1
    box = ones(nx, 1)*box;
end

if size(box, 2) == 1
    box_out = [-box box];
else
    assert(all(box(:, 2) >= box(:, 1)));
    box_out = box;
end

box_center = (box_out(:, 2) + box_out(:, 1))/2;
box_half   = (box_out(:, 2) - box_out(:, 1))/2;


end

