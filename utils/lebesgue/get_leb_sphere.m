%get a vector of lebesgue moments for a sphere
%in the YALMIP monomial basis
%partly based on getLebesgueMomentsNew (Milan Korda)
%author: Jared Miller, May 17, 2021
function moments = get_leb_sphere(d, n, r)
    %inputs:
    %   d:  maximum degree of moment sequence (0...d)
    %   n:  number of variables
    %   r:  radius of sphere    
    
    
    dv = monpowers(n, d);
    moments = LebesgueMomSphere(dv, r);
%     moments = zeros(size(dv,1),1);
%     r_scale = 1/(r^n);
%     for i = 1:numel(moments)
%         alpha = dv(i, :);
%         moments(i) = LebesgueMomSphere(alpha, r);
%     end
end