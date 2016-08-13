% A helper function that runs n4sid, calculates empirical eigenpairs

function [empvals, empvecs] = run_n4sid(data, noise, timestep, n, maxfreq, minfreq)

[len,num] = size(data);
range = max(max(data)) - min(min(data));
if(noise ~= 0)
    data = data + range*noise*rand(len,num)-1/2*range*noise;
    z = iddata(data,zeros(len,1),timestep);
    % set model order
    modelorder = n*2 + 4;
    m = n4sid(z, modelorder,'Form','modal','DisturbanceModel','estimate');
else
    z = iddata(data,zeros(len,1),timestep);
    % set model order
    modelorder = n*2 + 4;
    m = n4sid(z, modelorder,'Form','modal','DisturbanceModel','none');
end

%% Calculate Eigenvalue and Eigenvector Predictions from N4SID
[ mx, ~] = eig(m.A);
empvals = (log(eig(m.A))/timestep);
rangeactual = find(abs(imag(empvals)/2/pi) > minfreq & abs(imag(empvals)/2/pi) < maxfreq);
empvals = empvals(rangeactual);
[~, idx2] = sort(abs(imag(empvals)));
empvecs = m.C*mx;
empvecs = empvecs(:,rangeactual);
empvals = empvals(idx2);
empvecs = empvecs(:,idx2);

end