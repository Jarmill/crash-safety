figure(40)
clf
hold on
% fsurf(@(x, y) flow_func{4}.q([x; y]), [-1, 1, -1, 1]*box_lim);
SP = fsurf(@(x, y) subvalue_eval(x, y, flow_func), [-1, 1, -1, 1]*box_lim, 'MeshDensity', 150, 'EdgeColor', 'none');
% CP = fcontour(@(x, y) subvalue_eval(x, y, flow_func), [-1, 1, -1, 1]*box_lim, 'MeshDensity', 150);

%draw the unsafe set    %draw the moon
x_moon = moon_base(h_in, h_out);
moon_rot = [cos(moon_theta), sin(-moon_theta); sin(moon_theta), cos(moon_theta)];
Xu = moon_rot*x_moon*moon_scale + moon_center;
 

% Xu = Cu + circ_half* Ru;
patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')

C0 = [1; 0];
scatter3(C0(1), C0(2), 1, 100, 'k', 'filled')
c_c0 = subvalue_eval(C0(1), C0(2), flow_func)

xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
title('Flow Crash-Subvalue Map', 'fontsize', 14)
view(2)
axis square
cb = colorbar;
cb.Label.String='crash lower bound';
function c = subvalue_eval(x, y, flow_func)
    c = zeros(size(x));
    
    for i = 1:length(flow_func)
        c = max(c, flow_func{i}.q([x; y]));
    end
    c = min(c, ones(size(c)));
end