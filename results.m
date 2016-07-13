load 14bus/results.mat;
load 14bus/metadata.mat;
IMAGE = zeros(numcontigs);

for i = 1:numcontigs
   ranking = floor(tiedrank(M(i,:)));
   IMAGE(i,:) = ranking;
end

imagesc(IMAGE)