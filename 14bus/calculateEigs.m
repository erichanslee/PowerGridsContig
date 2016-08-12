%% Script takes data and calculates Eigs
%   More precisely, takes low rank updates and runs Arnoldi
%   with Sherman-Morrison-Woodbury calculations

%   For now uses eig, which should just be a QR iteration,
%   to calculate the eigenvalues.
load('metadata.mat');
A = matrix_read('data/matrix');

for i = 1:numcontigs
  n = length(A);
  u = matrix_read(sprintf('data/u%d', i));
  v = matrix_read(sprintf('data/v%d', i));
  [Q,H] = InvArnoldi(A,ones(n,1),n,0,u,v,sprank(u));
end
