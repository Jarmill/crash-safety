%Author: Milan Korda (I believe)
%originally titled momsphere.m

function y = LebesgueSphereMom(a,r)
% Returns the moment of monomial x^a on the unit hypersphere
% y = \int_{x : \sum_i x^2_i <= 1} x^a dx
% cf. Theorem 3.1 in J. B. Lasserre, E. S. Zeron. Solving a class
% of multivariate integration problems via Laplace techniques.
% Applicationes Mathematicae 28(4):391-405, 2001.


%INPUT:
% r - radius of the ball
% a - monomial powers


% The x1^a*x2^b lebesgue moment on a dim-dimensional sphere of
% radius r is r^{dim+a+b} * n-th moment on the unit sphere

% Example: Integrals of bi-variate monomials up to degree 10 over unit ball
% pows = monpowers(2,10);
% y = momsphere(pows,1)

if(~exist('r','var'))
    r = 1;
end

[m,n] = size(a);
y = zeros(m,1);
dim = size(a,2); % dimension

for k = 1:m

 if all(~rem(a(k,:),2))
  y(k) = prod(gamma((a(k,:)+1)/2))/ gamma(1+(n+sum(a(k,:)))/2);
  y(k) = r^(dim + sum(a(k,:))) * y(k);
 end

end


