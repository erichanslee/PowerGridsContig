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

function [actualvecs, empvecs, empresidual] = isolatecase(method,contignum, matrixnum, noise, window)

  maxfreq = .5;
  minfreq = .05;
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
  load('data/sim14_%d.mat', contignum);
  differential = DAE.n;
  algebraic = DAE.m;
  %A = matrix_read('matrix1');

  %%  Simulate Partial Placement of PMUs by obscuring percentage of simulated data

  rangebus = (DAE.n + Bus.n + 1):(DAE.n + Bus.n + Bus.n);
  PMUnum = round(window*length(rangebus));
  PMU = randsample(rangebus,PMUnum);
  PMU = sort(PMU);
  outrange = setdiff(rangebus, PMU);
  rangerest = [1:(DAE.n + Bus.n), (DAE.n + Bus.n + Bus.n + 1): (DAE.n + DAE.m), outrange];

  %% use n4sid
  data = Varout.vars(offset:end,PMU);
  m  = run_n4sid(data, noise, timestep, numlines);
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
  A = matrix_read(sprintf('data/matrixfull%d', contignum));
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

 % form DAE matrix E
I = eye(differential);
E = zeros(algebraic + differential);
E(1:differential,1:differential) = I;

for k = 1:numcontigs    
    %% Calculate Eigenvalue and Eigenvector Predictions from State Matrix
    % from the reduced state matrix
    A = full(matrix_read(sprintf('data/matrixfull%d', k)));
    format long
    
    %% Calculate Backward Error
    Ifull = eye(differential + algebraic);
    order = [PMU, rangerest];
    P = Ifull(order,:);
    out = zeros(length(temp2),1);
    for j = 1:length(temp2)

        % form variables to pass into calc_residual
        lambda = temp2(j);
        Ashift = A-lambda*E;
        xfull = zeros(differential + algebraic,1);
        x1 = actualvecs(:,j);
        res = calc_residual(method, Ashift, x1, PMU, rangerest, xfull, P);
        out(j) = norm(res);
    end
    data_dump(k) = mean(out);
end


end
