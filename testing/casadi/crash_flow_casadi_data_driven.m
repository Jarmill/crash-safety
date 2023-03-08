%optimal control of flow system using casadi
%try and get to the unsafe half-circle set
%use data-driven control inputs

%based on https://web.casadi.org/blog/ocp/ race_car.m

%crash bound 0.5499

N = 400; % number of control intervals

opti = casadi.Opti(); % Optimization problem

n = 2;
%% ---- data generation ------------
rng(33, 'twister')
Nsample = 40;
epsilon = [0; 0.5];
box_lim = 2;
Tmax = 5;
% epsilon = 2;
epsilon = [0; 0.5];
sample = struct('t', Tmax, 'x', @() box_lim*(2*rand(2,1)-1));
f_true= @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

DG = data_generator(sample);

x = sdpvar(2, 1);
mlist = monolist(x, 3);
model = struct('f0', [x(2);0], 'fw', [zeros(1, length(mlist)); mlist']);


observed = DG.corrupt_observations(Nsample, f_true, epsilon);


W = DG.data_cons(model, x, observed);
W.b = W.b - epsilon(2);
[w_handle, box]= DG.make_sampler(W);

L = size(W.A, 2);

f0 = polyval_func(model.f0, x);
fw = polyval_func(model.fw, x);

%% ---- decision variables ---------
X = opti.variable(n,N+1); % state trajectory
U = opti.variable(L,N);   % control trajectory (throttle)
T = opti.variable();      % final time
Z = opti.variable();

%% ---- dynamic constraints --------
% f = @(x,w) f_func(x) + [0; w]; % dx/dt = f(x,u)
f = @(x, w) f0(x) + fw(x)*w;

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
C0 = [1; 0];
opti.subject_to(X(1:2, 1) == C0);   % start at initial point
% opti.subject_to(X (3, 1) == Zmax);


% Cu = [1; -0.5];
Cu = [-0.25; -0.7];
Ru = 0.5;

%unsafe set 

c1f=@(x) Ru^2 - (x(1) - Cu(1)).^2 - (x(2) - Cu(2)).^2;

theta_c = 5*pi/4;
w_c = [cos(theta_c); sin(theta_c)];
c2f =@(x) w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 

Xu_con = @(x) [c1f(x); c2f(x)];
opti.subject_to(Xu_con(X(:, N+1)) >= 0);

%terminal time
opti.subject_to(5>=T>=0);
opti.subject_to(0<=Z<=Zmax);

%control effort
% opti.subject_to(-Z<=U<=Z);
for t = 1:size(U, 2)
    opti.subject_to(W.A*U(:, t)<=(Z+W.b));
end

%% ---- initial values for solver ---
opti.set_initial(T, 1);

%% ---- solve NLP              ------
opti.solver('ipopt'); % set numerical backend
sol = opti.solve();   % actual solve

fprintf('crash bound=%0.4f', sol.value(Z))

%% ---- plotting --------

figure(1)
clf
hold on

c = linspecer(3);

scatter(C0(1), C0(2), 200, 'ok')
plot(sol.value(X(1, :)), sol.value(X(2, :)), 'color', c(1, :), 'linewidth', 3)

theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
circ_half = [cos(theta_half_range); sin(theta_half_range)];
Xu = Cu + circ_half* Ru;
patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')

xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
title(sprintf('Crash States (z=%0.4f)', sol.value(Z)), 'fontsize', 16) 
pbaspect([diff(xlim), diff(ylim), 1])

figure(2)
clf
hold
% vtu = linspace(0, sol.value(T)-sol.value(dt), N);
vt = linspace(0, sol.value(T), N+1);
vtu = vt(1:end-1);
plot(vtu, sol.value(U))
% plot(xlim, sol.value(Z)*[1,1], 'k')
% plot(xlim, -sol.value(Z)*[1,1], 'k')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$u(t)$', 'interpreter', 'latex')
title('Data-Driven Crash Control-effort', 'Fontsize', 14) 

figure(3)
clf
hold
plot(vt, sol.value(X))
xlabel('$t$', 'interpreter', 'latex')
ylabel('$x(t)$', 'interpreter', 'latex')
title('Data-Driven Crash State', 'Fontsize', 14) 


figure(4)
clf
hold
plot(vtu, W.A*sol.value(U) - W.b)
plot(xlim, sol.value(Z)*[1,1], 'k')
plot(xlim, -sol.value(Z)*[1,1], 'k')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$\Gamma w(t)-h$', 'interpreter', 'latex')
title('Data-Driven Crash Control-Constraints', 'Fontsize', 14) 