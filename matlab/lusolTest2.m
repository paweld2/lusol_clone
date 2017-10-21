% lusolTest
% This is a script to test an implementation of the
% LUSOL Mex interface lu1fac.
% It loads data A and b suggested by Yin Zhang
% for debugging Kjell Eikland's C translation of LUSOL
% (with modifications by Yin Zhang).

% 14 Jan 2005: First version to test Yin's matrices.
%              They have m < n (but MINOS and SNOPT use only m>=n).
%              The trouble seems to be in the dense LU.
%              options.Dense2 = 1.1 disables the dense LU
%              and things then seem to work.

disp('TRP test')
%A = [1 1 3 2      % Yin's original example
%     1 1 1 1
%     0 0 1 0]

A = [1 1 1 1       % This one gives P=I, Q=I.
     1 3 2 1
     0 1 0 0]
[m,n] = size(A);

options            = lusolSet;
options.Pivoting   = 'TRP';
options.FactorTol  = 1.1;
 
% Gives L = [1          U = [1 1 1 1
%           [1 1               2 1 0
%             .5 1]              0 0]     Should be -.5 0]

options.Dense2 = 1.1;   % This disables switch to dense LU.
                        
[L,U,p,q,options] = lusolFactor(A,options);

inform = options.Inform;

if inform > 0
   disp(' ')
   disp('Hmmmm: luSOL(A) should have returned inform = 0.')
end

E = A - L*U;
e = norm(full(E),'fro');
disp(' ')
disp(['norm(A - L*U)_F = ' num2str(e)])

if e < 1e-8*(m*n)
  disp('This seems good')
else
  disp('This seems too large')
end
disp(' ')


disp('TCP test')
A = [1 1 2 1
     1 1 0 1
     5 5 0 1]
[m,n] = size(A);

options.Pivoting   = 'TCP';
options.FactorTol  = 1.1;
 
[L,U,p,q,options] = lusolFactor(A,options);

if inform > 0
   disp(' ')
   disp('Hmmmm: luSOL(A) should have returned inform = 0.')
end

E = A - L*U;
e = norm(full(E),'fro');
disp(' ')
disp(['norm(A - L*U)_F = ' num2str(e)])

if e < 1e-8*(m*n)
  disp('This seems good')
else
  disp('This seems too large')
end
disp(' ')
