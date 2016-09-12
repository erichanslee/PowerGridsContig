% Checks correctness of a few routines:

% place_PMU.m
% filter_eigenpairs.m
% id_contig.m
load metadata.mat

n = differential + algebraic;
A = rand(n);
E = eye(n);
method = 3;
[V,D] = eig(A,E);
win = 1:10;
empvecs = V(win,:);
empvals = diag(D);

[res,vecfull] = id_contig(A,E,method,empvals, empvecs,win);

if(norm(res) < 1e-8)
	disp('Function id_contig Passed\n');
else
	disp('Function id_contig Failed\n');
end

