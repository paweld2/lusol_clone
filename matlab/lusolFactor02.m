function [L,U,p,q,options] = lusolFactor(A,options)
%        options   = lusolSet;
%        [L,U,p,q,options] = lusolFactor(A,options);
%
% lusolFactor computes the sparse factorization A = L*U
% for a square or rectangular matrix A.  The vectors p, q
% are row and column permutations giving the pivot order.
%    L(p,p) is unit lower triangular with subdiagonals bounded
%           by options.FactorTol,
%    U(p,q) is upper triangular or upper trapezoidal,
%           depending on the shape of A.
%    A(p,q) = L(p,p)*U(p,q) would be a truly triangular
%             (or trapezoidal) permutation.
%    The rank of A tends to be reflected in the rank of U,
%    especially if options.Pivoting = 'TRP' or 'TCP'.
%
% Example use:
%    options = lusolSet(options,'pivot','TRP','factor',4.0);
% or equivalently
%    options = lusolSet;
%    options.Pivoting   = 'TRP';
%    options.FactorTol  = 4.0;

% The Matlab interface to LUSOL is maintained by
% Michael O'Sullivan and Michael Saunders, SOL, Stanford University.
%
% Known Bugs:
%
% The Fortran LUSOL is designed to any m-by-n A.
% The Mex interface lu1fac seems to work properly
% when m = n (square A).
%
% When m > n, L and U contain correct information SOMETIMES
% but they both have dimension m-by-n.
% lusolFactor could remove the last n-m rows of U to maintain A = L*U,
% but instead we add the last n-m columns of the identity to L.
%
% Beware, L*U is not always close to A for some reason.
% We keep the option m > n available in the hope that
% at least p and q are returned correctly.  They may be used
% to pull out a square submatrix of A, which can then be
% factored correctly.
%
% When m < n, the Mex interface returns an m-by-n L that
% crashes Matlab if it is accessed.  (L should be m-by-m.)
% lusolFactor issues a warning and returns.
% The permutations p and q may still be useful.
% Alternatively, consider factoring A' instead of A.
% With Rook Pivoting (options.Pivoting = 'TRP'), the same
% numerical properties are achieved.
%
% options.KeepLU = 'No' does not work.
% The Mex file still tries to create L and U.
%
% 02 Feb 1999: MJO: Developed MEX interface to LUSOL's lu1fac.
% 18 Oct 2000: MAS: Added options structure.
% 15 Apr 2001: MJO: options is now an output parameter.
% 14 Aug 2002: MAS: Added TRP option.
% 29 Apr 2004: MAS: L and U are returned as given by lu1fac.
%              There is no need to permute them to triangular form
%              because Matlab handles permuted triangles correctly
%              in statements like x = U\(L\b);

  
  [m,n] = size(A);
  if m < n
      disp(' ')
      disp('Warning: lu1fac is not safe with m < n.')
      disp('Refrain from accessing L and U.')
      disp('p and q may be ok.')
      disp('Also consider TRP on A(transpose).')
      disp(' ')
      % L = [];  U = [];  p = [];  q = [];
      % options.Inform = 100;     
  end

  if ~issparse(A), A = sparse(A); end     % Make sure A is sparse.

  luparm    = zeros(30,1);                % Store options in LUSOL arrays.
  parmlu    = zeros(30,1);

  luparm(1) = options.PrintFile;
  luparm(2) = options.PrintLevel;
  luparm(3) = options.MaxCol;
  if strcmpi(options.Pivoting,'TPP'), luparm(6) = 0; end
  if strcmpi(options.Pivoting,'TRP'), luparm(6) = 1; end
  if strcmpi(options.Pivoting,'TCP'), luparm(6) = 2; end
  if strcmpi(options.KeepLU  ,'No' ), luparm(8) = 0; end
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

  % 29 Apr 2004: No need to do this:
  % if options.KeepLU
  %        L = L(p,p);                  % Make L and U strictly triangular
  %        U = U(p,q);
  % end

  inform         = luparm(10);
  options.Inform = inform;
  options.Rank   = luparm(16);
  options.Nsing  = luparm(11);
  options.Growth = parmlu(16);

  if inform > 0
    disp(['Note: lu1fac returned Inform = ' num2str(inform)])
  end

  if m > n
    % U(p(n+1:m),:) = [];                 % Remove unwanted rows of U
    Im = speye(m);
    L  = [L Im(:,n+1:m)];
    disp(' ')
    disp('Note: m > n.  lu1fac may return a good p and q')
    disp('              but A = L*U may not hold')
  end
