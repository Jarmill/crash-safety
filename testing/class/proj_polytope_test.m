m = 5;
L = 6;
p = 2;

rng(4, 'twister')
A = randn(m, L);
G = randn(m, p);
b = randn(m, 1);

[A'; G']
