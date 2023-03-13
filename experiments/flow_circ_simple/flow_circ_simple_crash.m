%data driven crash-safety estimation of the flow system
%with a single input on x2dot

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
% epsilon = [0; 0.5];
sample = struct('t', Tmax, 'x', @() box_lim*(2*rand(2,1)-1));



%% generate model
yalmip('clear')
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
w_handle = @() wbound(2*rand()-1);
end
 
%% Solve SOS program
if SOLVE
    
    %start at a single point

%     C0 = [1.5; 0];
% %     C0 = [-1; 0];
%     R0 = 0.2;
C0 = [1; 0];
R0 = 0.4;

    INIT_POINT = 0;
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

%HAZARD (start in the unsafe set)
% X0 = Cu;


    theta_c = 5*pi/4;
    w_c = [cos(theta_c); sin(theta_c)];
    c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 

    Xu = struct('ineq', [c1f; c2f], 'eq', []);
    
    
    lsupp = loc_crash_options();
    lsupp.t = t;
    lsupp.TIME_INDEP = 0;
    lsupp.x = x;
    lsupp = lsupp.set_box(box_lim);
%     lsupp.X = struct('ineq', 2*box_lim^2 - sum(x.^2), 'eq', []);
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
    
    %Start from within an unsafe set;

    
    %INIT_POINT = 1
%  order=5; %    5.1179e-01
%    order=4; %5.0918e-01
 %  order=3;%4.3694e-01\
%  order=2; % 1.1843e-01
%     order=1;%1.1166e-07

%INIT_POINT = 0; casadi crash bound
order=1; %crash cost: 8.1010e-08
order=2; %crash cost: 6.5897e-02
% order=3; crash cost: 4.0539e-01
order=4; %crash cost: 4.6316e-01, time: 36.81
% order=5; %crash cost: 4.6381e-01, time: 868.81

    d = 2*order; 

    % [prog]= PM.make_program(d);
    % out = PM.solve_program(prog)
    out = PM.run(order);
    disp(sprintf('crash cost: %0.4e', out.obj))
    
    wbound = out.obj;
    w_handle = @() wbound*(2*rand()-1);
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

