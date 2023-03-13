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
        %unsafe set
    %moon heights
    h_in = 0.4;
h_out = 1;

    

%hugging the curve
moon_center = [0.4;-0.4];
moon_theta = -pi/10;
moon_scale = 0.8;

%same coordinates as half-circle example
% moon_center = [0;-0.7];
% moon_theta = -pi/4;
% moon_scale = 0.5

moon_rot = [cos(moon_theta), sin(-moon_theta); sin(moon_theta), cos(moon_theta)];
% x_moon_move = moon_rot*x_moon*moon_scale + moon_center;


%statistics of the moon
c_in = [0;0.5*(1/h_in - h_in)];
r_in = 0.5*(1/h_in + h_in);

c_out = [0;0.5*(1/h_out - h_out)];
r_out = 0.5*(1/h_out + h_out);

c_in_scale = moon_rot*c_in*moon_scale + moon_center;
c_out_scale = moon_rot*c_out*moon_scale + moon_center;

r_in_scale = moon_scale*r_in;
r_out_scale = moon_scale*r_out;

%constraints of the moon
con_inner =  sum((x-c_in_scale).^2) - r_in_scale^2;
con_outer =  -sum((x-c_out_scale).^2) + r_out_scale^2;
Xu = struct('ineq', [con_inner; con_outer], 'eq', []);
    
    
    %location support
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
%     lsupp.recover=1;
    lsupp.solver='mosek';
    lsupp.Zmax_Cap = 2;

    box = [-1, -1; 1, 1]*box_lim;
    lsupp.mom_handle = @(d) LebesgueBoxMom( d, box, 1);

    lsupp.verbose = 1;

%     objective = x(1);
   
    
    %% start up tester
    PM = crash_subvalue_sos(lsupp);

    
    %INIT_POINT = 1
    
    %order 5: primal infeasible. need to upper-bound the subvalue objective
    %to Zmax*vol(X)?


    %ZMax cap of 4:
%     order=1; %integral: 1.9777e-07, C0: 8.8079e-09
%     order=2;%integral: 1.4745e-07, C0: 6.8638e-10
%     order=3; %integral: 1.1393e+00, C0: -8.8452e-02
        order=4; %integral: 3.7723e+00, C0: -6.4243e-03

        %Zmax cap = 2
%     order=1; %integral: 1.3425e-07, C0: 4.6520e-10
%     order=2;%integral: 1.3425e-07, C0: 4.6520e-10
%     order=3; %integral: 1.0268e+00, C0: -7.8612e-02
%     order=4; %integral: 3.1875e+00, C0: -5.6924e-03
    order=5;

        
    d = 2*order; 

    
    out = PM.run(order);
    disp(sprintf('integral: %0.4e, C0: %0.4e', out.obj, out.func.q([1; 0])))
    
    load('subvalue_flow_moon_simple_2.mat', 'flow_func');
    flow_func{order} = out.func;
    save('subvalue_flow_moon_simple_2.mat', 'flow_func');
end

%% plot the subvalue function
if PLOT_SUBVALUE
    figure(40)
    clf
    hold on
    fsurf(@(x, y) flow_func{order}.q([x; y]), [-1, 1, -1, 1]*box_lim);
    
    %draw the moon
x_moon = moon_base(h_in, h_out);
moon_rot = [cos(moon_theta), sin(-moon_theta); sin(moon_theta), cos(moon_theta)];
Xu = moon_rot*x_moon*moon_scale + moon_center;
 

% Xu = Cu + circ_half* Ru;
patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')


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

