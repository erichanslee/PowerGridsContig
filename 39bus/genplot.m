function out = genplot()
load metadata.mat
M = zeros(numcontigs);
noise = 0;
window = 0;
for i = 1:numcontigs
    for j = 1:numcontigs
        %sanitycheck(i,j):
        % i = contingency number for simulation
        % j = matrix number i.e. state matrix
        M(i,:) = sanitycheck(i,noise,window);
    end
end

out = M;
end