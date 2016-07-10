%% Just a check to see whether or not everything matches up


clear all;
contignum = 3;
matrixnum = 3;
maxfreq = .5;
minfreq = .05;
initpsat;
load('metadata.mat')
%% run simulation, generate data
% at 60hz with fixed timesteps of 0.05s
%
Settings.freq = 60;
Settings.fixt = 1;
Settings.tstep = 0.05;

runpsat(strcat('contig',int2str(contignum)),'data');
runpsat('td');
for i = 1:numbuses
string = strcat('simulation/sim9bus' , num2str(i), '.txt');
fid = fopen(string, 'w');
fprintf(fid, '%f\n', Varout.vars(10:end,33+i));
end
%A = dlmread('matrix1'); A = spconvert(A); A = full(A);

%% use n4sid
data = Varout.vars(60:end,34:42);
[length,num] = size(data);
z = iddata(data,zeros(length,1),Settings.tstep);
m = n4sid(z, 20,'Form','modal','DisturbanceModel','none');

% n4sid gives discrete model with A_discrete = expm(A_cont*k)
% where k is the sampling time. 

[ mx, ~] = eig(m.A);

clc;
%% Calculate Eigenvalue and Eigenvector Predictions from State Matrix 
% from the reduced state matrix
I = eye(differential);
E = zeros(algebraic + differential);
E(1:differential,1:differential) = I;
predA = dlmread(strcat('data/matrixfull',int2str(matrixnum)));
predA = full(spconvert(predA));
[vi,di] = eig(predA,E); %solve generalized eigenvalue problem
predvalues = diag(di);

fprintf('Contingency %d simulated\n',contignum);
fprintf('Contingency %d predicted\n\n',matrixnum);


%% Sort and organize data properly
% actual = from system identification
% pred = from linearized system
temp1 = (diag(di));
temp2 = (log(eig(m.A))/Settings.tstep);

rangebus = (DAE.m + 1):(DAE.m + Bus.n);
rangepred = find(abs(imag(temp1)/2/pi) > .02 & abs(imag(temp1)/2/pi) <.4);  
rangeactual = find(abs(imag(temp2)/2/pi) > .02 & abs(imag(temp2)/2/pi) <.4);


temp1 = temp1(rangepred);
temp2 = temp2(rangeactual);

[~, idx1] = sort(abs(imag(temp1)));
[~, idx2] = sort(abs(imag(temp2)));

actualvecs = m.C*mx;
actualvecs = actualvecs(:,rangeactual);

%sort eigenvalues
temp1 = temp1(idx1);
temp2 = temp2(idx2);
%sort eigenvectors
predvecs = vi(rangebus,rangepred); %where the important eigenvectors are
predvecs = predvecs(:,idx1);
actualvecs = actualvecs(:,idx2);

format long
%% Check if eigenvectors match
fprintf('Checking if Eigenvectors Match up...\n');
fprintf('Column Index = Actual Eigenvectors\n');
fprintf('Row Index = Predicted Eigenvectors\n');
closeness = (normalizematrix(actualvecs)'*normalizematrix(predvecs));
printnorms(closeness)

%% Check if frequencies match
fprintf('Checking if Frequencies Match Up....\n');
fprintf('Predicted Frequencies from Linearized System\n');
fprintf('\n');
display(imag(temp1)/2/pi);
fprintf('Actual Frequencies from Simulation\n');
fprintf('\n');
display(imag(temp2)/2/pi);

%% Check if dampening matches
fprintf('Checking if Dampening Matches Up....\n');
fprintf('Predicted Frequencies from Linearized System\n');
fprintf('\n');
display(real(temp1));
fprintf('Actual Frequencies from Simulation\n');
fprintf('\n');
display(real(temp2));
%% Check if eigenvectors match

