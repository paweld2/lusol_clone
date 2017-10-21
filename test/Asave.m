function Asave(A)

%        Asave(A)
% saves the nonzeros of matrix A in file Afile.txt.
% It's probably better to delete Afile.txt first.

% 22 Jul 2007: First version of Asave.m, to save
%              Mike O'Sullivan's Q654.mat in a form
%              that lusoltest1.f can read.

  fid = fopen('Afile.txt','wt');
  [i,j,Aij] = find(A);
  [m,n] = size(A);
  nz    = length(Aij);
  
  for k=1:nz
    fprintf(fid,'%8i %8i %22.14e\n',i(k),j(k),Aij(k));
  end
  
  dir Afile.txt
  