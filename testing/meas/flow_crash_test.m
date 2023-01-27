if PLOT_DIST
    figure(2)
    clf
    hold on
    for i = 1:Nsample
        if i == 1
            plot(out_sim{i}.t, out_sim{i}.dist, 'c', 'DisplayName', 'Trajectories');
        else
            plot(out_sim{i}.t, out_sim{i}.dist, 'c', 'HandleVisibility', 'Off');
        end
    end
    
    plot(xlim, dist_rec*[1,1], '--r', 'LineWidth', 2 ,'DisplayName', 'Distance Bound');
    
    if optimal_pt
        plot(out_sim_peak.t, out_sim_peak.dist, 'b', 'LineWidth', 2, 'DisplayName','Closest Traj.');
        scatter(tp_rec, dist_rec, 300, '*b', 'LineWidth', 2,'DisplayName','Closest Point');
    end
    
    xlabel('time')
    ylabel('distance to unsafe set')
    title('Distance to unsafe set along trajectories', 'FontSize' , FS_title)
    legend('location', 'northwest')
end


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
