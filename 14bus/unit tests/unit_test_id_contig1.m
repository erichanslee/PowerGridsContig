% Checks correctness of id_contig.m

cd ..;

%% Test Case: Random (Works Fine)
load metadata.mat
maxfreq = 10;
minfreq = .5;
n = differential + algebraic;
A = rand(n);
E = rand(n);
method = 3;
[v1,d1] = eig(A,E);
[v1, d1] = filter_eigpairs(minfreq, maxfreq, diag(d1), v1);

win = (differential + numlines + 1):(differential + numlines + numlines);
v1_win = v1(win,:);
d1_diag = diag(d1);
v1_arg = normalizematrix(v1_win);

[res1,vec1] = id_contig(A,E,method,d1_diag, v1_arg,  win);

disp('Dot Product of Fitted and Real eigenvector:');
disp(abs(vec1(:,1)'*v1(:,1)/norm(v1(:,1))));

cd 'unit tests'/