function w = rej_sample_poly(A, b, box)
%REJ_SAMPLE_POLY Rejection sampling to generate uniformly-distributed 
%data inside a convex polytope A w <= b
%
%go to hit and run later

[m, n] = size(A);

if nargin <= 2
    box = ones(n, 1)*[-1, 1];
end

done = 0;
iter = 0;
while ~done
    w = rand(n, 1).*(box(:, 2)-box(:, 1)) + box(:, 1);
    
    if all(A*w <= b)
        done = 1;
    end
    iter = iter + 1;
end


end

