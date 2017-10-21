function [L, U, p, q] = splu(A)

options = lusolSet;
options.Pivoting = 'TCP';

[L, U, p, q] = lusolFactor(A, options);
