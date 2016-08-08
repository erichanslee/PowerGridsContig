%% Checks whether or not contingency matches reduced state matrix modes
% contignum = the contingency to simulate
% noise = white noise amount (input as percentage the max amplitude)
% window = percentage of full view into power network
%           implemented by randomly removing subset of voltage readings
function [predcontig, actualcontig] = testbed(noise, window)



maxfreq = .5;
minfreq = .05;
initpsat;
load('metadata.mat')
contignum = ceil(numcontigs*rand);
actualcontig = contignum;
%% Basic Pre-Run Checks
if(noise > 1 || noise < 0)
    error('Problems with parameter "Noise". Please enter in a real number in the range of [0,1]')
end

if(window > 1 || window < 0)
    error('Problems with parameter "Window". Please enter in a real number in the range of [0,1]')
end

if(contignum > numcontigs)
    error('Contigency Number not found!')
end

%% run simulation, generate data
% at 60hz with fixed timesteps of 0.05s, simulating 20 polls/second on PMU
%
Settings.freq = 60;
Settings.fixt = 1;
Settings.tstep = 0.05;
offset = 50;

runpsat(strcat('contig',int2str(contignum)),'data');
runpsat('td');
for i = 1:numbuses
    string = strcat('simulation/sim9bus' , num2str(i), '.txt');
    fid = fopen(string, 'w');
    fprintf(fid, '%f\n', Varout.vars(offset:end,DAE.n + Bus.n + i));
    fclose('all');
end

differential = DAE.n;
algebraic = DAE.m;
%A = dlmread('matrix1'); A = spconvert(A); A = full(A);

%%  Simulate Partial Placement of PMUs by obscuring percentage of simulated data

rangebus = (DAE.n + Bus.n + 1):(DAE.n + Bus.n + Bus.n);
PMUnum = round(window*length(rangebus));
PMU = randsample(rangebus,PMUnum);
PMU = sort(PMU);
outrange = setdiff(rangebus, PMU);
rangerest = [1:(DAE.n + Bus.n), (DAE.n + Bus.n + Bus.n + 1): (DAE.n + DAE.m), outrange];
N = length(PMU);

%% use n4sid
data = Varout.vars(offset:end,PMU);
size(data)
[len,num] = size(data);
range = max(max(data)) - min(min(data));
if(noise ~= 0)
    data = data + range*noise*rand(len,num)-1/2*range*noise;
    z = iddata(data,zeros(len,1),Settings.tstep);
    % set model order
    modelorder = Bus.n*2 + 4;
    m = n4sid(z, modelorder,'Form','modal','DisturbanceModel','estimate');
else
    z = iddata(data,zeros(len,1),Settings.tstep);
    % set model order
    modelorder = Bus.n*2 + 4;
    m = n4sid(z, modelorder,'Form','modal','DisturbanceModel','none');
end

% n4sid gives discrete model with A_discrete = expm(A_cont*k)
% where k is the sampling time.

%%  Print some more stuff
clc;

%% Calculate Eigenvalue and Eigenvector Predictions from N4SID

[ mx, ~] = eig(m.A);
temp2 = (log(eig(m.A))/Settings.tstep);
rangeactual = find(abs(imag(temp2)/2/pi) > minfreq & abs(imag(temp2)/2/pi) < maxfreq);
temp2 = temp2(rangeactual);
[~, idx2] = sort(abs(imag(temp2)));
actualvecs = m.C*mx;
actualvecs = actualvecs(:,rangeactual);
temp2 = temp2(idx2);
actualvecs = actualvecs(:,idx2);
actualvecs = normalizematrix(actualvecs);

data_dump = zeros(1,numcontigs);
%% Calculate Eigenvalue and Eigenvector Predictions from State Matrix
% from the reduced state matrix
I = eye(differential);
E = zeros(algebraic + differential);
E(1:differential,1:differential) = I;
A = dlmread(strcat('data/matrixfull',int2str(contignum)));
A = full(spconvert(A));
[vi,di] = eig(A,E); %solve generalized eigenvalue problem

%% Sort and organize data from State Matrix properly
% actual = from system identification
% pred = from linearized system
temp1 = (diag(di));
rangepred = find(abs(imag(temp1)/2/pi) > minfreq & abs(imag(temp1)/2/pi) < maxfreq);
temp1 = temp1(rangepred);
[~, idx1] = sort(abs(imag(temp1)));

%sort eigenvalues
temp1 = temp1(idx1);

%sort eigenvectors
predvecs = vi(rangebus,rangepred); %where the important eigenvectors are
predvecs = predvecs(:,idx1);
predvecsEntire = vi(:,rangepred);
predvecsEntire = predvecsEntire(:,idx1);

for k = 1:numcontigs
    
    %% Calculate Eigenvalue and Eigenvector Predictions from State Matrix
    % from the reduced state matrix
    I = eye(differential);
    E = zeros(algebraic + differential);
    E(1:differential,1:differential) = I;
    A = dlmread(strcat('data/matrixfull',int2str(k)));
    A = full(spconvert(A));
    
    
    
    format long
    
    %% Calculate Backward Error
    out = zeros(length(temp2),1);
    Ifull = eye(DAE.n + DAE.m);
    order = [PMU, rangerest];
    P = Ifull(order,:);
    for j = 1:length(temp2)
        
        % Form the shifted matrix
        lambda = temp2(j);
        Ashift = A-lambda*E;
        
        % Solve an OLS problem to fill in unknown entries (min residual)
        xfull = zeros(DAE.n + DAE.m, 1);
        xfull(PMU) = actualvecs(:,j);
        xfull(rangerest) = Ashift(:,rangerest)\(Ashift(:,PMU)*xfull(PMU));
        
        % Compute the residual and save the norm
        res = Ashift*xfull;
        out(j) = norm(res);
        
    end
    data_dump(k) = mean(out);
    
end
[~, idx] = min(data_dump);
predcontig = idx;


end
