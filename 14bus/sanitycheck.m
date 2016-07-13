%% Just a check to see whether or not everything matches up
function out = sanitycheck(contignum, matrixnum)

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
offset = 50;

runpsat(strcat('contig',int2str(contignum)),'data');
runpsat('td');

differential = DAE.n;   
algebraic = DAE.m;
%A = dlmread('matrix1'); A = spconvert(A); A = full(A);
rangebus = (DAE.n + Bus.n + 1):(DAE.n + Bus.n + Bus.n);

%% use n4sid
data = Varout.vars(offset:end,rangebus);
[length,num] = size(data);
z = iddata(data,zeros(length,1),Settings.tstep);
% set model order
modelorder = Bus.n*2 + 4;
m = n4sid(z, modelorder,'Form','modal','DisturbanceModel','none');

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

%organizing data
rangepred = find(abs(imag(temp1)/2/pi) > minfreq & abs(imag(temp1)/2/pi) < maxfreq);  
rangeactual = find(abs(imag(temp2)/2/pi) > minfreq & abs(imag(temp2)/2/pi) < maxfreq);


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
closeness = printnorms(closeness);


%% Check if frequencies match
fprintf('Checking if Frequencies Match Up....\n');
fprintf('Predicted Frequencies from Linearized System\n');
fprintf('\n');
disp(imag(temp1)/2/pi);
fprintf('Actual Frequencies from Simulation\n');
fprintf('\n');
disp(imag(temp2)/2/pi);

%% Check if dampening matches
fprintf('Checking if Dampening Matches Up....\n');
fprintf('Predicted Frequencies from Linearized System\n');
fprintf('\n');
disp(real(temp1));
fprintf('Actual Frequencies from Simulation\n');
fprintf('\n');
disp(real(temp2));

%% Using data to calculate measure of closeness
[rowsize, colsize] = size(closeness);
sum = 0;
if(colsize > rowsize)
    for i = 1:rowsize
        [value, idx] = max(closeness(i,:));
        sum = sum + value;
        closeness(:,idx) = 0;
    end
elseif (colsize < rowsize)
    for i = 1:colsize
        [value, idx] = max(closeness(:,i));
        sum = sum + value;
        closeness(idx,:) = 0;
    end
elseif (colsize == rowsize)
    for i = 1:colsize
        [value, idx] = max(closeness(:,i));
        sum = sum + value;
        closeness(idx,:) = 0;
    end
else
    error('Problem with Inner Product Size');
end


out = sum;
end