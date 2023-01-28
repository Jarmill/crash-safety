m = 10;
L = 6;
p = 0;

rng(4, 'twister')
A = randn(m, L);
G = randn(m, p);
b = ones(m, 1);

P = [A'; G']


N = null(P);

x = sdpvar(size(N, 2), 1);

sol = optimize([norm(x, 2)<= 1; N*x>=0], []);

x_rec = value(x);
Nx_rec = value(N*x);