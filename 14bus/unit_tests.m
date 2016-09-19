% Checks correctness of id_contig.m

%% Test Case: Random (Works Fine)
load metadata.mat

% n = differential + algebraic;
% A = rand(n);
% E = eye(n);
% method = 3;
% [v1,d1] = eig(A,E);
% win = (differential + numlines + 1):(differential + numlines + numlines);
% v1_win = v1(win,:);
% d1_diag = diag(d1);
% 
% [res1,vec1] = id_contig(A,E,method,d1_diag, normalizematrix(v1_win), win);


%% Test Case: Power System (Doesn't work well)

maxfreq = .5;
minfreq = .05;
contignum = 1;
n = differential + algebraic;
method = 3;

I = eye(differential);
E = zeros(algebraic + differential);
E(1:differential,1:differential) = I;
A = full(matrix_read(sprintf('data/matrixfull%d', contignum)));
win = (differential + numlines + 1):(differential + numlines + numlines);
[v2,d2] = eig(A,E); 
[v2_subset, d2_subset] = filter_eigpairs(minfreq, maxfreq, diag(d2), v2);

% Pick out window and normalize %
v2_arg = normalizematrix(v2_subset(win,:));
[res2, vec2] = id_contig(A, E, method, d2_subset, v2_arg, win);

%% 
disp('Dot Product of Fitted and Real eigenvector in Second Problem:');
disp(vec2(:,2)'*v2(:,1)/norm(v2(:,1)));