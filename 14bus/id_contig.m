% ~~~~~~~~~INPUTS~~~~~~~~~ %

% method = method type number one would like to use
% A,E = generalized eigenvalue problem matrices
% empvals = empirical eigenvalues
% empvecs = empirical eigenvectors
% PMU = indices x1 is located on

% ~~~~~~~~~OUTPUTS~~~~~~~~~ %
% fittedvecs = fitted eigenvectors 

function [fittedres, fittedvecs] = id_contig(A, E, method, empvals, empvecs, PMU)
    load metadata.mat

    fittedvecs = zeros(differential + algebraic, length(empvals));
    for j = 1:length(empvals)

        % form variables to pass into calc_residual
        lambda = empvals(j);
        Ashift = A-lambda*E;
        xfull = zeros(differential + algebraic,1);
        x1 = empvecs(:,j);
        rangerest = 1:(differential + algebraic);
        rangerest = rangerest(~ismember(rangerest, PMU));
        [fittedres(:,j), fittedvecs(:,j)] = calc_residual(method, Ashift, x1, PMU, rangerest, xfull);
    end
end

% ~~~~~~~~~INPUTS~~~~~~~~~ %

% method = method type number one would like to use
% Ashift = matrix in question
% x1 = subset of eigenvector
% PMU = indices x1 is located on
% rangerest = indices of the rest of the eigenvector
% xfull = empty vector of size length(PMU) + length(rangerest)
% P = permutation matrix passed in for  

% ~~~~~~~~~OUTPUTS~~~~~~~~~ %

% residual = calculated residual
% vec = full fitted eigenvector

function [residual, vec] = calc_residual(method, Ashift, x1, PMU, rangerest, xfull)
        switch method
            case 1	%% METHOD 1
                % Solve an OLS problem to fill in unknown entries (min residual)
                xfull(PMU) = x1;
                xfull(rangerest) = (-1*Ashift(:,rangerest))\(Ashift(:,PMU)*xfull(PMU));
                
                % Compute the residual and save the norm
                residual = Ashift*xfull;
                
                
            case 2	%% METHOD 2
                % Solve an OLS problem to fill in unknown entries (min residual)
                xfull(PMU) = x1;
                xfull(rangerest) = (-1*Ashift(:,rangerest))\(Ashift(:,PMU)*xfull(PMU));
                
                % Compute the residual and save the norm
                xfull = xfull/norm(xfull);
                residual= Ashift*xfull;
                
                
            case 3  %% METHOD 3
                Ifull = eye(differential + algebraic);
                order = [PMU, rangerest];
                P = Ifull(order,:);

                Ashift = Ashift*ctranspose(P);
                
                % Form Gramian
                T = zeros(DAE.n + DAE.m,1+length(rangerest));
                T(1:length(PMU),1) = x1;
                T((length(PMU)+1):end,2:end) = eye(length(rangerest));
                G = ctranspose(T)*(ctranspose(Ashift)*Ashift)*T;
                
                % Calculate smallest eigenvector and then form eigenvector
                [vs,ds] = eigs(G,1,'sm');
                xfull(1:length(PMU)) = vs(1)*x1;
                xfull((length(PMU)+1):end) = vs(2:end);
                
                % Compute the residual and save the norm
                
                residual = Ashift*xfull;
                % Note: could easily just use eigenvalue as output
                % but we want the full eigenvector for debugging purposes
                
            case 4  %% METHOD 4: Making x1 unit length again
                Ifull = eye(differential + algebraic);
                order = [PMU, rangerest];
                P = Ifull(order,:);
                Ashift = Ashift*ctranspose(P);
                % Form Gramian
                T = zeros(DAE.n + DAE.m,1+length(rangerest));
                T(1:length(PMU),1) = x1;
                T((length(PMU)+1):end,2:end) = eye(length(rangerest));
                G = ctranspose(T)*(ctranspose(Ashift)*Ashift)*T;
                
                % Calculate smallest eigenvector and then form eigenvector
                [vs,ds] = eigs(G,1,'sm');
                xfull = zeros(DAE.n + DAE.m,1);
                xfull(1:length(PMU)) = vs(1)*x1;
                xfull((length(PMU)+1):end) = vs(2:end);
                
                % Compute the residual, renormalize x1 to have unit length
                % and save the norm
                residual = 1/vs(1)*Ashift*xfull;
        end
    vec = xfull; 
end