% Checks correctness of id_contig.m

%% Test Case: Random
load metadata.mat

n = differential + algebraic;
A = rand(n);
E = eye(n);
method = 3;
[v1,d1] = eig(A,E);
win =  (differential + numlines + 1):(differential + numlines + numlines);
v1_win = v1(win,:);
d1_diag = diag(d1);

[res1,vec1] = id_contig(A,E,method,d1_diag, normalizematrix(v1_win), win);


%% Test Case: Power System

maxfreq = .5;
minfreq = .05;
contignum = 1;

I = eye(differential);
E = zeros(algebraic + differential);
E(1:differential,1:differential) = I;
A = full(matrix_read(sprintf('data/matrixfull%d', contignum)));
[v2,d2] = eig(A,E); 
d2 = d2 + 1e-6*eye(size(E));


[v2_subset, d2_subset] = filter_eigpairs(minfreq, maxfreq, diag(d2), v2);

[res2, vec2] = id_contig(A, E, method, d2_subset, normalizematrix(v2_subset(win,:)), win);