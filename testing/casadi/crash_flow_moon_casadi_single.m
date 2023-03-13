%optimal control of flow system using casadi
%try and get to the unsafe half-circle set
%use an Linf control penalty

%based on https://web.casadi.org/blog/ocp/ race_car.m

N = 400; % number of control intervals

opti = casadi.Opti(); % Optimization problem

n = 2;

%% ---- decision variables ---------
X = opti.variable(n,N+1); % state trajectory
U = opti.variable(1,N);   % control trajectory (throttle)
T = opti.variable();      % final time
Z = opti.variable();

%% ---- dynamic constraints --------
f_func = @(x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
f = @(x,w) f_func(x) + [0; w]; % dx/dt = f(x,u)

dt = T/N; % length of a control interval
for k=1:N % loop over control intervals
   % Runge-Kutta 4 integration
   k1 = f(X(:,k),         U(:,k));
   k2 = f(X(:,k)+dt/2*k1, U(:,k));
   k3 = f(X(:,k)+dt/2*k2, U(:,k));
   k4 = f(X(:,k)+dt*k3,   U(:,k));
   x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4); 
   opti.subject_to(X(:,k+1)==x_next); % close the gaps
end

%% ---- objective          ---------

% Z = X(3, end);
% J = sum(U.^2)*dt;
opti.minimize(Z);
% opti.minimize(T); % race in minimal time

%% ---- boundary conditions --------


Zmax = 1;
Tmax = 5;
%initial point
% C0 = [1.5; 1]; %cool 0.
% C0 = [1.25; 1];
% C0 = [-1; 0]; 
C0 = [0; 0]; %0.3232
% C0 = [1; 0]; %0.2240
% C0 = [0.5; 0];
% C0 = [0.4;-0.4]; %the coolest, has a twist, which means its a local min
opti.subject_to(X(1:2, 1) == C0);   % start at initial point
% opti.subject_to(X (3, 1) == Zmax);

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
con_inner =  @(x) sum((x-c_in_scale).^2) - r_in_scale^2;
con_outer = @(x) -sum((x-c_out_scale).^2) + r_out_scale^2;
%unsafe set 
Xu_con = @(x) [con_inner(x); con_outer(x)];
opti.subject_to(Xu_con(X(:, N+1)) >= 0);

%terminal time
opti.subject_to(5>=T>=0);
opti.subject_to(0<=Z<=Zmax);

%control effort
opti.subject_to(-Z<=U<=Z);

%% ---- initial values for solver ---
opti.set_initial(T, 1);

%% ---- solve NLP              ------
opti.solver('ipopt'); % set numerical backend
sol = opti.solve();   % actual solve

fprintf('crash bound=%0.4f\n', sol.value(Z))

%% ---- plotting --------

figure(1)
clf
hold on

% scatter(C0(1), C0(2), 200, 'ok')
scatter(sol.value(X(1, 1)), sol.value(X(2, 1)), 100, c(1, :), 'filled')
c = linspecer(3)
plot(sol.value(X(1, :)), sol.value(X(2, :)), 'color', c(1, :), 'linewidth', 2)
% plot(sol.value(X(1, :)), sol.value(X(2, :)))

%draw the moon
x_moon = moon_base(h_in, h_out);
moon_rot = [cos(moon_theta), sin(-moon_theta); sin(moon_theta), cos(moon_theta)];
Xu = moon_rot*x_moon*moon_scale + moon_center;
 

% Xu = Cu + circ_half* Ru;
patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')

xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
title(sprintf('Crash States (z=%0.4f)', sol.value(Z)), 'fontsize', 16) 
pbaspect([diff(xlim), diff(ylim), 1])

figure(2)
clf
hold
vt = linspace(0, sol.value(T), N);
plot(vt, sol.value(U))
plot(xlim, sol.value(Z)*[1,1], 'k')
plot(xlim, -sol.value(Z)*[1,1], 'k')
xlabel('t')
ylabel('u(t)')
title('Crash Control-effort') 
