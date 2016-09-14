function out = normalizematrix(A)
%% function just takes a matrix and normalizes each column. 

[~,width] = size(A);
for i = 1:width
   A(:,i) = A(:,i)/norm(A(:,i)); 
end

out = A;
end
