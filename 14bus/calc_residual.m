%	function taking in matrix and subset of potential eigenvector
%	and outputting a residual subject to method number

function [res, vec] = calc_residual(method, Ashift, x1, PMU, rangerest, xfull, P)
        switch method
            case 1	%% METHOD 1
                
                % Solve an OLS problem to fill in unknown entries (min residual)
                xfull(PMU) = x1;
                xfull(rangerest) = (-1*Ashift(:,rangerest))\(Ashift(:,PMU)*xfull(PMU));
                
                % Compute the residual and save the norm
                res = Ashift*xfull;
                
                
            case 2	%% METHOD 2
                
                % Solve an OLS problem to fill in unknown entries (min residual)
                xfull(PMU) = x1;
                xfull(rangerest) = (-1*Ashift(:,rangerest))\(Ashift(:,PMU)*xfull(PMU));
                
                % Compute the residual and save the norm
                xfull = xfull/norm(xfull);
                res = Ashift*xfull;
                
                
            case 3  %% METHOD 3
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
                
                res = Ashift*xfull;
                % Note: could easily just use eigenvalue as output
                % but we want the full eigenvector for debugging purposes
                
            case 4  %% METHOD 4: Making x1 unit lengh again
                
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
                res = 1/vs(1)*Ashift*xfull;
        end
    vec = xfull; 
end