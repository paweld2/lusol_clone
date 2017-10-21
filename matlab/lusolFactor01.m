function [L,U,p,q,options] = lusolFactor(A,options)
%        options   = lusolSet
%        [L,U,p,q,options] = lusolFactor(A,options);
%
% lusolFactor computes the sparse factorization A(p,q) = L*U,
% where  p and q are row and column permutations,
%        L is unit lower triangular with subdiagonals bounded
%          by options.FactorTol,
%        U is upper triangular (or upper trapezoidal).
% The rank of A tends to be reflected in the rank of U.
%
% Example change of parameters:
%        options = lusolSet(options,'pivot','TCP','factor',5.0);
% or     options.Pivoting  = 'TCP';
%        options.FactorTol = 5.0;

% The Matlab interface to LUSOL is maintained by
% Michael O'Sullivan and Michael Saunders, SOL, Stanford University.
%
% 02 Feb 1999: MJO: Developed MEX interface to LUSOL's lu1fac.
% 18 Oct 2000: MAS: Added options structure.
% 15 Apr 2001: MJO: Added output options, options now output parameter
% 14 Aug 2002: MAS: Added TRP option.

  if ~issparse(A), A = sparse(A); end     % Make sure A is sparse.

  luparm    = zeros(30,1);                % Store options in LUSOL arrays.
  parmlu    = zeros(30,1);

  luparm(1) = options.PrintFile;
  luparm(2) = options.PrintLevel;
  luparm(3) = options.MaxCol;
  if strcmpi(options.Pivoting,'TRP'), luparm(6) = 1; end
  if strcmpi(options.Pivoting,'TCP'), luparm(6) = 2; end
  if strcmpi(options.KeepLU  ,'Yes'), luparm(8) = 1; end

  parmlu(1) = options.FactorTol;
  parmlu(2) = options.UpdateTol;
  parmlu(3) = options.DropTol;
  parmlu(4) = options.Utol1;
  parmlu(5) = options.Utol2;
  parmlu(6) = options.Uspace;
  parmlu(7) = options.Dense1;
  parmlu(8) = options.Dense2;

  [L,U,p,q,luparm,parmlu] = lu1fac(A,luparm,parmlu);   % Factorize A = L*U

  if options.KeepLU
          L = L(p,p);                    % Make L and U strictly triangular
          U = U(p,q);
  end

  [m, n] = size(A);  
  options.Inform = luparm(10);
  options.Rank   = luparm(16);
  options.Nsing  = luparm(11);
  options.Growth = parmlu(16);
