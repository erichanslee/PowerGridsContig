%   orthprojection projects a vector x along the column space of a matrix
%   A via a QR and then a projection.

%   INPUT:
%   x = the vector you want to project
%   A = the matrix whose column space you want to use
%   mode = boolean checking whether A is already orthgonal
%   out = projection residual

function out = orthprojection(x,A,mode)
  if(mode == 1)
    Q = A;
    [~,len] = size(Q);
    proj = zeros(size(x));
    for i = 1:len
      q = Q(:,i);
      proj = proj + dot(x,q)*q;
    end
  else
    Q = orth(A);
    [~,len] = size(Q);
    proj = zeros(size(x));
    for i = 1:len
      q = Q(:,i);
      proj = proj + dot(q,x)*q;
    end
  end
  out = x - proj;
end
