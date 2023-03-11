%data driven crash-safety estimation of the flow system
%with a single input on x2dot

%break up the sections here into functions

PROBLEM = 1;
SOLVE = 1;
SAMPLE = 0;
PLOT = 0;
PLOT_SUBVALUE = 1;

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
% epsilon = [0; 0.5];
sample = struct('t', Tmax, 'x', @() box_lim*(2*rand(2,1)-1));



%% generate model
t = sdpvar(1, 1);
x = sdpvar(2, 1);

model = struct('f0', f_true(0, x), 'fw', [0; 1]);

W = struct('A', [1; -1], 'b', [0; 0], 'G', []);

% W = DG.data_cons(model, x, observed);
% W_true = W;
% W.b = W.b - epsilon(2);
% W.G = [];

%TODO: write code to prune constraints/possibly trivial faces

%modify W for the crash program
% [model_cheb,W_cheb] = DG.center_cheb(model, W);
% W_red = DG.reduce_constraints(W_cheb);
% [w_handle, box]= DG.make_sampler(W_true);

w_handle = @() 2*rand()-1;
end
 
%% Solve SOS program
if SOLVE
    
    %start at a single point

    
    %unsafe set
    
theta_c = 3*pi/2;

Cu = [1; -0.5];
Ru = 0.5;
    c1f = Ru^2 - sum((x-Cu).^2);
%     c2f = -diff(x-Cu);

    
    w_c = [cos(theta_c); sin(theta_c)];
    c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 

    Xu = struct('ineq', [c1f; c2f], 'eq', []);
    
    
    lsupp = loc_crash_options();
    lsupp.t = t;
    lsupp.TIME_INDEP = 0;
    lsupp.x = x;
    lsupp = lsupp.set_box(box_lim);
    lsupp.X_init = lsupp.X;
    lsupp.X_term = Xu;
    lsupp.f0 = model.f0;
    lsupp.fw = model.fw;
    lsupp.Tmax = Tmax;
    lsupp.W = W;
    lsupp.recover=1;
    lsupp.Zmax = 2;
    lsupp.solver='mosek';

    box = [-1, -1; 1, 1]*box_lim;
    lsupp.mom_handle = @(d) LebesgueBoxMom( d, box, 1);

    lsupp.verbose = 1;

%     objective = x(1);

    objective = [c1f; c2f];
    
    %% start up tester
    PM = crash_subvalue_sos(lsupp);

    
    %INIT_POINT = 1
    

%     order=1; %integral 1.0424e-07,  X01 eval: 3.9574e-09, X02 eval: 6.5930e-09 
%     order =2;  integral: 3.1026e-07, X01 eval: -6.4925e-09, X02 eval: -5.2559e-09 
%     order=3; %subvalue integral: 1.7682e+00, X01 eval: -1.0750e-01, X02 eval: -5.2519e-02 
%     order=4; %subvalue integral: 4.2754e+00,  X01 eval: -8.6549e-02, X02 eval: -6.9769e-02 
    order=5;
    


    d = 2*order; 

    
    out = PM.run(order);
    X01 = [0.395; 1.21];
    X02 = [1.279; -1.21];

    fprintf('subvalue integral: %0.4e \n X01 eval: %0.4e \n X02 eval: %0.4e \n', ...
        out.obj, out.func.q(X01), out.func.q(X02))
    
    load('subvalue_same_dist.mat', 'flow_func');
    flow_func{order} = out.func;
    save('subvalue_same_dist.mat', 'flow_func');
end

%% plot the subvalue function
if PLOT_SUBVALUE
    figure(40)
    clf
    hold on
    fsurf(@(x, y) flow_func{order}.q([x; y]), [-1, 1, -1, 1]*box_lim);
    
    %draw the unsafe set
    theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
    circ_half = [cos(theta_half_range); sin(theta_half_range)];
    Xu = Cu + circ_half* Ru;
    patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none')
    
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

