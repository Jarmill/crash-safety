function [x_moon] = moon_base(h_in, h_out, Narc)
%MOON_BASE set of points on a moon
%The moon contains points [-1, 0] and [1, 0]. It is defined between a
%circle at height h_in and at height h_out, where the height is the lowest
%y coordinate on the circle.


if nargin < 3
    Narc = 200;
end

%centers and radii of the circle
c_in = 0.5*(1/h_in - h_in);
r_in = 0.5*(1/h_in + h_in);

c_out = 0.5*(1/h_out - h_out);
r_out = 0.5*(1/h_out + h_out);

if h_in == 0
    x_moon_in = [1,-1; 0, 0];
else
    angle_moon_in  = asin(1/r_in);
    if c_in >= 0
        theta_moon_in  = linspace(angle_moon_in, -angle_moon_in, Narc)-pi/2;
    else
        theta_moon_in  = linspace(2*pi-angle_moon_in, angle_moon_in, Narc)+pi/2;
    end
    x_moon_in = (r_in)*[cos(theta_moon_in); sin(theta_moon_in)] + [0; c_in];
end

angle_moon_out = asin(1/r_out);
if c_out >=0
    theta_moon_out = linspace(-angle_moon_out, angle_moon_out, Narc)-pi/2;
else
    theta_moon_out  = linspace(angle_moon_out, 2*pi-angle_moon_out, Narc)+pi/2;
end

x_moon_out = (r_out)*[cos(theta_moon_out); sin(theta_moon_out)] + [0; c_out];

% x_in_refl = x_in
x_moon = [x_moon_out, x_moon_in];
% x_moon = x_moon_out;

end

