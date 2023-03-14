f_func = @(x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
Tmax = 5;
BOX = 2;

%initial set
C0 = [1.5; 0];
R0 = 0.4;




%% draw the vector field

Ngrid = 21;
xx = linspace(-BOX, BOX, Ngrid);
[XX, YY] = meshgrid(xx);

UU = zeros(Ngrid);
VV = zeros(Ngrid);
for i = 1:Ngrid
    for j = 1:Ngrid
        fcurr = f_func([xx(j), xx(i)]);
        UU(i, j) = fcurr(1);
        VV(i, j) = fcurr(2);
    end
end

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


%% sample a pair of trajectories

% X0 = [0.395, 1.21; ...
%     1.279, -1.21];

% X0 = [0.3948, 1.21; ...
%     1.279, -1.21];

% X0 = [0, 1; 1.297 -1.5];
% X0 = [0, 1; 1.2966 -1.5];
C01 = [0, 1];
C02 = [1.2966 -1.5];

X0 = {C01, C02};


NP = size(X0, 2);

supp_curr = @(t, x) box_event(t, x, BOX);

options = odeset('Events',supp_curr, 'RelTol', 1e-9);
osm = cell(NP, 1);
dist_close = zeros(NP, 1);
for i = 1:NP
    osm{i} = ode23(@(t, x) f_func(x), [0, Tmax], X0{i}, options);
    dist_close_curr = zeros(size(osm{i}.x));
    for k = 1:length(dist_close_curr)
        dist_close_curr(k) = aff_half_circ_dist(osm{i}.y(:, k), Ru, theta_c, Cu);        
    end
    osm{i}.dist_close_vec = dist_close_curr;
    dist_close(i)= min(dist_close_curr);
end

%% plot the field
Ntheta = 100;
x_dist = dist_contour(Ntheta, Ru, max(dist_close)) + Cu;

cc = linspecer(4);

figure(2)
clf
% tiledlayout(1, 2)

% t1 = nexttile;
hold on
streamslice(XX, YY, UU, VV)

patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')
axis equal
xlim([-BOX, BOX])
ylim([-BOX, BOX])

for i = 1:NP
    plot(osm{i}.y(1, :), osm{i}.y(2, :), 'LineWidth', 3, 'color', cc(2+i, :));
    X0c = X0{i};
    scatter(X0c(1), X0c(2), 100, cc(2+i, :), 'filled')
end
% scatter(X0(:, 1), X0(:, 2), 100, 'k')
plot(x_dist(1, :), x_dist(2, :), 'LineWidth', 3, 'color', 'r');

xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')

% t2 = nexttile;
% hold on
% quiver(XX, YY, UU, VV, 4)
% 
% Ntheta = 100;
% x_dist = dist_contour(Ntheta, Ru, max(dist_close)) + Cu;
% % streamline(stream2(XX, YY, UU, VV, XX, YY))
% % streamslice(XX, YY, UU, VV)
% 
% patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')
%    axis equal
% 
% for i = 1:NP
%     plot(osm{i}.y(1, :), osm{i}.y(2, :), 'LineWidth', 3, 'color', cc(2+i, :));
% end
% 
% plot(x_dist(1, :), x_dist(2, :), 'LineWidth', 3, 'color', 'r');
% 
% scatter(X0(:, 1), X0(:, 2), 100, 'k')
% 
% xlim([-BOX, BOX])
% ylim([-BOX, BOX])
% 

%% function helpers
function dist_out = half_circ_dist(x_in, R)
    %return the L2 distance  between the point x_in and the half circle
    %||x_in||^2 <= R^2 intersect x_in(2) <= 0.
%     reshape(x_in, [], 1);
    if x_in(2) >= 0
        %flat region
        if x_in(1) < -R
            dist_out = hypot(x_in(1)+R, x_in(2));
        elseif x_in(1) > R
            dist_out = hypot(x_in(1)-R, x_in(2));
        else
            dist_out = x_in(2);
        end
    else
        %circle region
        dist_out = max(norm(x_in, 2)-R, 0);
    end

end

function dist_out = aff_half_circ_dist(x_in, R, theta_c, Cu)

    theta_cf = theta_c - 3*pi/2;
    Rot_mat = [cos(theta_cf) -sin(theta_cf); sin(theta_cf) cos(theta_cf)];
    x_aff = Rot_mat'*(x_in - Cu);
    
    dist_out = half_circ_dist(x_aff, R);

end
   
function x_dist = dist_contour(Ntheta, R, c)
    %compute a contour at distance c away from the half-circle with N_theta
    %sample points


    theta_q1 = linspace(0, pi/2, Ntheta);
    theta_q2 = linspace(pi/2, pi, Ntheta);
    theta_q34 = linspace(pi, 2*pi, 2*Ntheta);

    %contour level
    

    x_right = [c*cos(theta_q1)+R; c*sin(theta_q1)];
    x_left = [c*cos(theta_q2)-R; c*sin(theta_q2)];
    x_bottom = [(c+R)*cos(theta_q34); (c+R)*sin(theta_q34)];
    x_dist = [x_right, x_left, x_bottom];
end
