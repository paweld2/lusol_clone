function [L,U,P,Q] = lu1fac(A,pivot,tolerance,memscalar)
%
% [L,U,P,Q] = LU1FAC(A,pivot,tol,memscalar) is functionally similar to
% [L,U,P,Q] = lu(A).  It is implemented as a wrapper around the LUSOL
% library.
%
% file:      	lu1fac.m
% directory:    /u/yzhang/MATLAB/mxLUSOL/
% created: 	Fri Nov 25 2005 
% author:  	Yin Zhang 
% email:   	yzhang@cs.utexas.edu
%

  if nargin < 2, pivot = 'trp';     end
  if nargin < 3, tolerance = 5;     end
  if nargin < 4, memscalar = 2;     end
  
  % it must match thos defined in lusol.h
  switch(lower(pivot))
   case {'tpp'}
    pivotmethod = 0;
   case {'trp'}
    pivotmethod = 1;
   case {'tcp'}
    pivotmethod = 2;
   case {'tsp'}
    pivotmethod = 3;
   otherwise
    pivotmethod = 1; % trp
  end

  if (~issparse(A))
    A = sparse(A);
  end
  
  [m,n] = size(A);
  Im = speye(m,m);
  In = speye(n,n);
  
  % XXX: lusol/mxlu1fac used to fail if there are empty cols in the beginning
  % or the middle of A.  To alleviate the problem, we sort A based on sum(A~=0)
  [ignore,q0] = sort(-sum(A~=0,1));
  A  = A(:,q0);
  [ignore,p0] = sort(-sum(A~=0,2));
  A  = A(p0,:);
  
  [L,U,p,q] = mxlu1fac(A,pivotmethod,tolerance,memscalar);
  
  % compute P
  P  = Im(p0(p),:);

  % compute Q (need to include the effect of corder)
  Q = In(:,q0(q));
  
  % make L lower triagular
  L  = Im - L(p,p);
  
  % make U upper triangular
  U  = U(q,:)';


