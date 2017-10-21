function [T1] = lutridi( n, d )
%        [T1] = lutridi( n, d );
% Illustrates instability in LU.
%
  
  e  = ones(n,1);
  T  = spdiags([e d*e e], (-1:1), n, n);
  T1 = [ T(2:n,:)      % Permute row 1 to bottom
	 T(1  ,:) ];
  [L,U] = lu(T1,0);    % Force diagonal pivots
  b  = T1*e;
  x  = U\(L\b);
  error = norm(x-e)
  keyboard
  
  