% genscript(basename, numbuses, numcontigs)
%
% Generate a directory of basic power flow data for all contingencies.
%
% Inputs:
%     basename: Base name (no extension) of PSAT spec file
%     numbuses: Number of buses in PSAT spec
%     numcontigs: Number of contingencies to consider
%
% Output files:
%     data/matrixfull:  Full Jacobian for base problem
%     data/matrixfullN: Full Jacobian for contingency N
%     data/matrixred:   Reduced Jacobian for base problem
%     data/matrixredN:  Reduced Jacobian for contingency N
%     data/uN:          Left factor for base-contingency reduced J
%     data/vN:          Right factor for "
%
function genscript(basename, numbuses, numcontigs)

  % Compute full and redued Jacobians for base case and write to file
  % DSB: Why do we trim tiny elements here and not below?
  [A, Aoriginal] = jacobian(basename);
  Aoriginal = trim_tiny(Aoriginal);
  write_sparse('data/matrixfull', A);
  write_sparse('data/matrixred', Aoriginal);

  for i = 1:numcontigs

    % Set up data file for contingency
    remove_line('contig.m', sprintf('%s.m', basename), i+1);

    % Compute full and reduced Jacobians and write to file
    [A, As] = jacobian('contig');
    write_sparse(sprintf('data/matrixfull%d', i), A);
    write_sparse(sprintf('data/matrixred%d', i), As);

    % Compute SVD and write factors to file
    [U, S, V] = svds(Aoriginal-As);
    write_sparse(sprintf('data/u%d', i), trim_tiny(U*S));
    write_sparse(sprintf('data/v%d', i), trim_tiny(V));

    % Clean up data file
    delete('contig.m');

  end
end


% ------------------------------------------------------
% Remove an indicated line from file infname and write to outfname
%
function remove_line(outfname, infname, linenum)
  str = fileread(infname);
  pos = strfind(str, sprintf('\n'));
  str(pos(linenum-1)+1 : pos(linenum)) = [];
  fid = fopen(outfname, 'w');
  fprintf(fid, '%s', str);
  fclose(fid);
end


% ------------------------------------------------------
% Remove small entries from a matrix
%
function [A] = trim_tiny(A, thresh)
  if nargin < 2, thresh = 1e-10; end
  A(abs(A) < thresh) = 0;
end


% ------------------------------------------------------
% Run PSAT power flow to get full and reduced Jacobian matrices
%
% DSB: Why are we subtracting 1e-6*I?  For that matter, why are
%      we explicitly forming the Schur complement at all?
%
function [A, Ared] = jacobian(fname)
  initpsat;
  runpsat(fname, 'data');
  runpsat('pf');
  A = [DAE.Fx DAE.Fy; DAE.Gx DAE.Gy];
  Ared = DAE.Fx - DAE.Fy*(DAE.Gy\DAE.Gx)- 1e-6*speye(DAE.n);
end


% ------------------------------------------------------
% Write a sparse matrix to a file
function write_sparse(fname, A)
  [col,row,val] = find(A);
  data_dump = [col,row,val];
  dlmwrite(fname, data_dump, 'precision', 16);
end
