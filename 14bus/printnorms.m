%%prints the norms of a matrix entry wise

function out = printnorms(A)
    [n,m] = size(A);
    N = zeros(n,m);
    for i = 1:n
        for j = 1:m
            N(i,j) = norm(A(i,j));
        end
    end
    
    disp(N);
    out = N;
end