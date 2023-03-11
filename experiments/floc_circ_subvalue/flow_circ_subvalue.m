%data driven peak estimation of a linear system

%break up the sections here into functions

PROBLEM = 1;
SOLVE = 1;
SAMPLE = 0;
PLOT = 1;

if PROBLEM
rng(33, 'twister')
%% generate samples
% A_true = [-1 1; -1 -0.3];
f_true = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
n=2;

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
    lsupp.X_init = lsupp.X;
    lsupp.X_term = Xu;
    lsupp.f0 = model.f0;
    lsupp.fw = model.fw;
    lsupp.Tmax = Tmax;
    lsupp.W = W;    
    lsupp.solver='mosek';
    
%     box_supp = box_process(2, box_lim);
%     lsupp.mom_handle = @(d) LebesgueBoxMom(d, box_supp', 1);
    lsupp.mom_handle = @(d) get_leb_sphere(d, n, sqrt(2)*box_lim);
    
    lsupp.verbose = 1;

%     objective = x(1);

    objective = [c1f; c2f];
    
    %% start up tester
    PM = crash_subvalue_sos(lsupp);

    
    flow_func = cell(4, 1);
    
    %INIT_POINT = 1
    
%int q(x) dx = 11.6027, almost entirely constant q=0.461655113926 (WAS A
%BUG!)
%true cost at C0:  0.54999
%     order = 1; %integral: 2.1926e-01, C0: 0.006180037245630
%     order = 2; %integral: 3.8185e+00, C0:  0.1829
    order=3; %integral: 7.8326e+00, C0: 0.3399
%   order=4; %integral: , C0: 

    d = 2*order; 

    % [prog]= PM.make_program(d);
    % out = PM.solve_program(prog)
    out = PM.run(order);
    disp(sprintf('crash integral: %0.4e', out.obj))
    
    
    load('subvalue_flow_circ.mat', 'flow_func');
    flow_func{order} = out.func;
    save('subvalue_flow_circ.mat', 'flow_func');
    
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

    sample_time = toc;
end

%% plot trajectories
if PLOT
    
%     if PLOT

%at order 3: @(varnew)[0.467354323477+0.000167921538097.*varnew(1,:)+0.00169771773626.*varnew(2,:)-0.000498411376915.*varnew(1,:).^2+0.00133310886064.*varnew(1,:).*varnew(2,:)+0.000539686485551.*varnew(2,:).^2+0.000798228685011.*varnew(1,:).^3+0.000198782740172.*varnew(1,:).^2.*varnew(2,:)+0.000121361325098.*varnew(1,:).*varnew(2,:).^2+1.00263239178e-06.*varnew(2,:).^3+0.000745764745296.*varnew(1,:).^4+0.000926683580285.*varnew(1,:).^3.*varnew(2,:)+0.000110698242649.*varnew(1,:).^2.*varnew(2,:).^2+0.000128896943176.*varnew(1,:).*varnew(2,:).^3-8.26819841086e-06.*varnew(2,:).^4-9.05489261143e-05.*varnew(1,:).^5-3.39839771908e-07.*varnew(1,:).^4.*varnew(2,:)+2.50365260024e-06.*varnew(1,:).^3.*varnew(2,:).^2-3.42647650485e-06.*varnew(1,:).^2.*varnew(2,:).^3-1.45159899011e-06.*varnew(1,:).*varnew(2,:).^4+7.85818593753e-09.*varnew(2,:).^5+3.58353128595e-06.*varnew(1,:).^6-1.07105951638e-06.*varnew(1,:).^5.*varnew(2,:)-1.63871259057e-06.*varnew(1,:).^4.*varnew(2,:).^2-4.0870780829e-06.*varnew(1,:).^3.*varnew(2,:).^3-1.80255411768e-06.*varnew(1,:).^2.*varnew(2,:).^4-1.80622167422e-06.*varnew(1,:).*varnew(2,:).^5+1.31309703981e-06.*varnew(2,:).^6]
figure(55)
clf
fsurf(@(x, y) out.func.q([x; y]), [-box_lim, box_lim, -box_lim, box_lim]);
    
    
%     PS = peak_sos_plotter(out, out_sim);

%     DG.data_plot_2(observed);
%     PS.state_plot_2(box_lim);
%     if INIT_POINT
%         scatter(C0(1), C0(2), 200, 'k')
%     else
%         theta = linspace(0,2*pi, 200);
%         plot(R0*cos(theta)+C0(1), R0*sin(theta)+C0(2), 'color', 'k', 'LineWidth', 3);
%     end
%     
%     %plot the unsafe set
%     theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
%     circ_half = [cos(theta_half_range); sin(theta_half_range)];
%     Xu = Cu + circ_half* Ru;
%     patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none')
%     
    %observation plot    


    
    
end

