% testing

%path(path, '..\..\..\version3.6')
%path(path, '..\..\..\factorize')

%load Q_654

thresh = 0.5;
drop = eps;
utol1 = eps^(2/3);
utol2 = eps^(2/3);

options = lusolSet;
options.Pivoting = 'TCP';
options.DropTol = drop;
options.FactorTol = 1 / thresh;
options.PrintLevel = 0;
options.Utol1 = utol1;
options.Utol2 = utol2;

[L, U, p, q, options] = lusolFactor(Q_c, options);
