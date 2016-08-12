% [out] = normalizematrix(A)
%
% Scale matrix columns to unit two-norm.
%
function out = normalizematrix(A)
  out = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
end
