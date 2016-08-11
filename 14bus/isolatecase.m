% isolatecase is an auxillary function for debugging purposes that checks
% runs a time-domain simulation for contingency <contignum> and then
% calculates residuals for the linearized system <matrixnum>

% That is, it checks data from contingency identification routine on a
% cases-by-case basis for debugging purposes.

% ~~~~~~~~~INPUTS~~~~~~~~~ %

% method = method type number one would like to use
% contignum = contingency to simulate
% matrixnum = linearized system to use for residual calculation
% noise = percentage of max amplitude to add as gaussian noise
% window = percentage of PMUs visible

function [proportion] = isolatecase(method,contignum, matrixnum, noise, window)



maxfreq = .5;
minfreq = .05;
initpsat;
load('metadata.mat')

%% Basic Pre-Run Checks
if(method > 5 || method < 1)
    error('Problems with parameter "Method". Please an integer in [1,3]')
end

if(noise > 1 || noise < 0)
    error('Problems with parameter "Noise". Please enter in a real number in the range of [0,1]')
end

if(window > 1 || window < 0)
    error('Problems with parameter "Window". Please enter in a real number in the range of [0,1]')
end

if(contignum > numcontigs)
    error('Contigency Number not found! Please enter a smaller integer')
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

%% Sort and organize data from State Matrix properly (NOT NEEDED FOR NOW)
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

k = matrixnum;

%% Calculate Eigenvalue and Eigenvector Predictions from State Matrix
% from the reduced state matrix
I = eye(differential);
E = zeros(algebraic + differential);
E(1:differential,1:differential) = I;
A = dlmread(strcat('data/matrixfull',int2str(k)));
A = full(spconvert(A));



format long

%% Calculate Backward Error
Ifull = eye(DAE.n + DAE.m);
order = [PMU, rangerest];
P = Ifull(order,:);
out = zeros(length(temp2),1);
proportion = zeros(1,length(temp2));
if( length(temp2) ~= length(temp1) && method == 5)
    return
end
for j = 1:length(temp2)
    switch method
        case 1	%% METHOD 1
            % Form the shifted matrix
            lambda = temp2(j);
            Ashift = A-lambda*E;
            
            % Solve an OLS problem to fill in unknown entries (min residual)
            xfull1 = zeros(DAE.n + DAE.m, 1);
            xfull1(PMU) = actualvecs(:,j);
            xfull1(rangerest) = (-1*Ashift(:,rangerest))\(Ashift(:,PMU)*xfull1(PMU));
            
            % Compute the residual and save the norm
            res = Ashift*xfull1;
            out(j) = norm(res);
            proportion(j) = norm(xfull1(PMU))/norm(xfull1);
            
        case 2	%% METHOD 2
            lambda = temp2(j);
            Ashift = A-lambda*E;
            
            xfull2 = zeros(DAE.n + DAE.m, 1);
            xfull2(PMU) = actualvecs(:,j);
            xfull2(rangerest) = (-1*Ashift(:,rangerest))\(Ashift(:,PMU)*xfull2(PMU));
            
            % Normalize as well
            xfull2 = xfull2/norm(xfull2);
            res = Ashift*xfull2;
            out(j) = norm(res);
            proportion(j) = norm(xfull2(PMU));
            
        case 3  %% METHOD 3
            % Form the shifted matrix
            lambda = temp2(j);
            Ashift = (A-lambda*E)*P';
            
            % Form Gramian
            T = zeros(DAE.n + DAE.m,1+length(rangerest));
            T(1:length(PMU),1) = actualvecs(:,j);
            T((length(PMU)+1):end,2:end) = eye(length(rangerest));
            G = T'*(Ashift'*Ashift)*T;
            
            % Calculate smallest eigenvector and then form eigenvector
            [vs,ds] = eigs(G,1,'sm');
            xfull3 = zeros(DAE.n + DAE.m,1);
            xfull3(1:length(PMU)) = vs(1)*actualvecs(:,j);
            xfull3((length(PMU)+1):end) = vs(2:end);
            
            % Compute the residual and save the norm
            res = Ashift*xfull3;
            out(j) = norm(res);
            proportion(j) = abs(vs(1));
            % ~~~Note~~~: could easily just use eigenvalue as output
            % but we want the full eigenvector for debugging purposes
            
        case 4  %% METHOD 4: Making x1 unit lengh again
            % Form the shifted matrix
            lambda = temp2(j);
            Ashift = (A-lambda*E)*P';
            
            % Form Gramian
            T = zeros(DAE.n + DAE.m,1+length(rangerest));
            T(1:length(PMU),1) = actualvecs(:,j);
            T((length(PMU)+1):end,2:end) = eye(length(rangerest));
            G = T'*(Ashift'*Ashift)*T;
            
            % Calculate smallest eigenvector and then form eigenvector
            [vs,ds] = eigs(G,1,'sm');
            xfull3 = zeros(DAE.n + DAE.m,1);
            xfull3(1:length(PMU)) = vs(1)*actualvecs(:,j);
            xfull3((length(PMU)+1):end) = vs(2:end);
            
            % Compute the residual, renormalize x1 to have unit length
            % and save the norm
            res = 1/vs(1)*Ashift*xfull3;
            out(j) = norm(res);
            proportion(j) = abs(vs(1));
            
        case 5 %% Not a Method, simply checking theoretical eigenvectors
            
            proportion(j) = norm(predvecs(:,j))/norm(predvecsEntire(:,j));
    end
    
    
end

display('done');





end
