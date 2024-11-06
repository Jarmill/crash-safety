function [observed] = corrupt_observations(Np,sample, f, epsilon)
%Generate Nsample corrupted observations of the system x'(t) = f(t,x)
%with an L_infinity bounded noise of epsilon
%Inputs:
%   Np:         number of sampled points to run
%   sampler:    struct with fields 't' and 'x' of function_handles
%   f:          system dynamics
%   epsilon:    noise bound
%
%Output:        struct 'observed'
%                   t:          time points
%                   x:          state points
%                   xdot_true:  noiseless observations
%                   xdot_noise: noisy observations
%                   epsilon:    noise bound


%% initial fill of output


x0 = pull_x();
nx = length(x0);
observed = struct('t', zeros(1, Np), 'x', zeros(nx, Np), ...
    'xdot_true', zeros(nx, Np), 'xdot_noise', zeros(nx, Np), 'epsilon', epsilon);


noise_sample = (2*rand(nx, Np)-1)*epsilon;


%% generate sampled points

if isnumeric(sample.t)
    Tmax = sample.t;
    sample.t = @() rand()*Tmax;
end


for i = 1:Np
   tcurr = sample.t();
   xcurr = pull_x();
   
   
   xdotcurr = f(tcurr, xcurr);
   
   
   observed.t(i) = tcurr;
   observed.x(:, i) = xcurr;
   observed.xdot_true(:, i) = xdotcurr;
%    observed.xdot_noise(:, i) = xdotcurr + noise_sample(:, i);
end

observed.xdot_noise = observed.xdot_true + noise_sample;

    function x0 = pull_x()
        if isa(sample.x, 'function_handle')
            x0 = sample.x();      
        else
            %numeric
            %assume that the dimensions of x and w are compatible if arrays are
            %passed in.
            if size(sample.x, 2) == 1
                %single point
                x0 = sample.x;
            else
                %array, so Ns = size(sampler, 2)
                x0 = sample.x(:, i);
            end
        end
    end


end

