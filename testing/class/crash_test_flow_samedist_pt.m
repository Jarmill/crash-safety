%support sets

t = sdpvar(1, 1);
x = sdpvar(2, 1);

X0_pt =[0.3948, 1.21;
    1.279, -1.21];

%point (0.3948, 1.21): order4 = 0.3769
%point (1.279, -1.21): order4 = 0.6028


%unsafe set
theta_c = 3*pi/2;
w_c = [cos(theta_c); sin(theta_c)];
Cu = [1; -0.5];
Ru = 0.5;
c1f = Ru^2 - (x(1) - Cu(1)).^2 - (x(2) - Cu(2)).^2;
c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 

Xu = struct('ineq', [c1f; c2f], 'eq', 0);

INIT_POINT = 1;
box_lim =  [-0.6, 1.75;-1.5,1.5];
Tmax = 5;
Zmax = 2;





f_true = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

model = struct;
model.f0 = f_true(t, x);
model.fw = [0; 1];
poly = struct('A', [1; -1], 'b', [1; 1], 'G', []);

%% options
lsupp = loc_crash_options();
lsupp.t = t;
lsupp.TIME_INDEP = 0;
lsupp.x = x;
lsupp = lsupp.set_box(box_lim);
lsupp.X_term = Xu;
lsupp.f0 = model.f0;
lsupp.fw = model.fw;
lsupp.poly = poly;
lsupp.Tmax = Tmax;
lsupp.Zmax = Zmax;
lsupp.verbose = 1;

%% model setup

lsupp.X_init = X0_pt(2, :)';

order = 4;
PM = crash_sos(lsupp);
out = PM.run(order);
disp(out.obj)
% lsupp.get_all_supp()
% XT = lsupp.get_X_term()