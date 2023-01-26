mset clear
mpol('w', 1, 1);
mpol('z', 1, 1);
mpol('x', 2, 1);
mpol('t', 1, 1);

vars = struct;
vars.t = t;
vars.x = [x; z];
vars.w = w;


f_func = @(x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
Tmax = 5;
BOX = 2;
ZMAX = 2;



f = [f_func(x) + [0; w]; 0];

%initial set
C0 = [1.5; 0];
R0 = 0.4;


%% sample a pair of trajectories

% X0 = [0.4, 1.2; 1.25, -1.2];
X0 = [0.395, 1.21; ...
    1.279, -1.21];

%% location support 

lsupp = loc_support(vars);
lsupp = lsupp.set_box([-BOX, BOX; -BOX, BOX; 0, ZMAX]);
% lsupp.X_init = X0;
lsupp.Tmax = 5;
lsupp.FREE_TERM = 0;
lsupp.disturb = [w<=z; -w<=z];
objective = z;

%% describe the unsafe set
%unsafe set
% theta_c = 5*pi/4;      %safe, bound = 0.3184
theta_c = 3*pi/2;

Cu = [1; -0.5];
Ru = 0.5;
%unsafe set 
theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
circ_half = [cos(theta_half_range); sin(theta_half_range)];
Xu = Cu + circ_half* Ru;

c1f = Ru^2 - (x(1) - Cu(1)).^2 - (x(2) - Cu(2)).^2;

% theta_c = 3*pi/2;
% theta_c = 
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 

zcon = z*(ZMAX - z)>=0;

lsupp.X_term = [[c1f; c2f] >= 0;zcon];


%% solve the program 
order_list = 1:3; 
NP = size(X0, 2);
sol = cell(NP, length(order_list));
z_opt = zeros(NP, length(order_list));
status= zeros(NP, length(order_list));
for i = 1:NP
    for k = 1:length(order_list)
        lsupp.X_init = [x(1:2)==[X0(i, 1); X0(i, 2)]; zcon];
        PM = peak_manager(lsupp, f, objective);
        PM.MAXIMIZE = 0;
        
        sol{i, k} = PM.run(order_list(k), Tmax);
        z_opt(i, k) = sol{i, k}.obj_rec;
        status(i, k) = sol{i, k}.status;
    end
end


   
