%data driven peak estimation of a linear system

%break up the sections here into functions

PROBLEM = 1;
SOLVE = 1;
SAMPLE = 0;
PLOT = 0;

if PROBLEM
rng(33, 'twister')
%% generate samples
% A_true = [-1 1; -1 -0.3];
f_true = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];


% Nsample = 100;
% Nsample = 50;
Nsample = 40;
% Nsample = 30;
% Nsample = 20;
% Nsample = 4;
box_lim = 2;
Tmax = 5;
% epsilon = 2;
epsilon = [0; 0.5];
sample = struct('t', Tmax, 'x', @() box_lim*(2*rand(2,1)-1));



%% generate model
t = sdpvar(1, 1);
x = sdpvar(2, 1);

DG = data_generator(sample);

observed = DG.corrupt_observations(Nsample, f_true, epsilon);
% [model, W] = DG.reduced_model(observed, x, 1, 1);
% model = DG.poly_model(vars, 3);
mlist = monolist(x, 3);
model = struct('f0', [x(2);0], 'fw', [zeros(1, length(mlist)); mlist']);

W = DG.data_cons(model, x, observed);
W_true = W;
W.b = W.b - epsilon(2);
W.G = [];

%TODO: write code to prune constraints/possibly trivial faces

%modify W for the crash program
% [model_cheb,W_cheb] = DG.center_cheb(model, W);
% W_red = DG.reduce_constraints(W_cheb);
[w_handle, box]= DG.make_sampler(W_true);

end
 
%% Solve SOS program
if SOLVE
    
    %start at a single point

%     C0 = [1.5; 0];
% %     C0 = [-1; 0];
    R0 = 0.2;
C0 = [1; 0];
% R0 = 0.4;
    INIT_POINT = 1;
    if INIT_POINT
        X0 = C0;
    else
        X0 = struct('ineq', R0^2 - sum((x-C0).^2), 'eq', []);
    end
    
    %unsafe set
    
%     Ru = 0.3;
%     Cu = [0; -0.5];
Cu = [-0.25; -0.7];
Ru = 0.5;
    c1f = Ru^2 - sum((x-Cu).^2);
%     c2f = -diff(x-Cu);

    theta_c = 5*pi/4;
    w_c = [cos(theta_c); sin(theta_c)];
    c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 

    Xu = struct('ineq', [c1f; c2f], 'eq', []);
    
    
    lsupp = loc_crash_options();
    lsupp.t = t;
    lsupp.TIME_INDEP = 0;
    lsupp.x = x;
    lsupp = lsupp.set_box(box_lim);
    lsupp.X = struct('ineq', 2*box_lim^2 - sum(x.^2), 'eq', []);
    % lsupp = lsupp.set_box(3);
    lsupp.X_init = X0;
    lsupp.X_term = Xu;
    lsupp.f0 = model.f0;
    lsupp.fw = model.fw;
    lsupp.Tmax = Tmax;
    lsupp.W = W;
    lsupp.recover=0;
    lsupp.solver='mosek';

    
    lsupp.recover=0;
    
    lsupp.verbose = 1;

%     objective = x(1);

    objective = [c1f; c2f];
    
    %% start up tester
    PM = crash_sos(lsupp);

    
    %INIT_POINT = 1
    
    %BAD BOUNDS:
%     order = 4;% 6.2094e-01, taking 4979.66 seconds = 1.38 hours
    %but this violates the crash-bound from crash_flow_casadi_data_driven.
    %why? Bad conditioning? Invalid solution? BECAUSE MY CODE WAS BUGGED!
    %Gram0 has a condition number of 4.5327e+08 (if that matters)
%     order=3; %  5.0660e-01 (slight infeasibility, gram eig = 4*10^-5. doesn't work on sdpa_gmp, too big
%     order = 2; %  4.6207e-01
%     order = 1; %4.6166e-01

    %GOOD BOUNDS
%     order=4; %5.4999e-01 (99.38 minutes)
order=3; %4.8649e-01 something. need to re-run
    %     order=2; %4.4231e-01
%     order=1; %5.8215e-02


    d = 2*order; 

    % [prog]= PM.make_program(d);
    % out = PM.solve_program(prog)
    out = PM.run(order);
    disp(sprintf('crash cost: %0.4e', out.obj))
    
end

%% Sample trajectories
if SAMPLE
    
    s_opt = sampler_options;
    s_opt.mu = 0.2;
    if INIT_POINT
        s_opt.sample.x = @() X0;
    else
        s_opt.sample.x = @() R0*ball_sample(1,2)'+C0;
    end
    s_opt.sample.d = w_handle;
    s_opt.Nd = size(model.fw, 2);
    
    s_opt.Tmax = lsupp.Tmax;
    s_opt.parallel = 1;
    
%     Nsample_traj = 10;
    Nsample_traj = 120;
    
    tic
    out_sim = sampler(out.dynamics, Nsample_traj, s_opt);

%     out_sim = traj_eval(out, out_sim);

    sample_time = toc;
end

%% plot trajectories
if PLOT
    
%     if PLOT
    
    PS = peak_sos_plotter(out, out_sim);
%     PS.nonneg_zeta();
%     PS.obj_plot();
%     PS.state_plot();
%     PS.v_plot();
%     PS.nonneg_traj();
    DG.data_plot_2(observed);
    PS.state_plot_2(box_lim);
    if INIT_POINT
        scatter(C0(1), C0(2), 200, 'k')
    else
        theta = linspace(0,2*pi, 200);
        plot(R0*cos(theta)+C0(1), R0*sin(theta)+C0(2), 'color', 'k', 'LineWidth', 3);
    end
    
    %plot the unsafe set
    theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
    circ_half = [cos(theta_half_range); sin(theta_half_range)];
    Xu = Cu + circ_half* Ru;
    patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none')
    
    %observation plot    


    
    
end

