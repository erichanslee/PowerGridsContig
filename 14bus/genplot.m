function out = genplot()
load metadata.mat
M = zeros(numcontigs);

for i = 1:numcontigs
    for j = 1:numcontigs
        %sanitycheck(i,j):
        % i = contingency number for simulation
        % j = matrix number i.e. state matrix
        M(i,j) = sanitycheck(i,j);
    end
end

out = M;
end