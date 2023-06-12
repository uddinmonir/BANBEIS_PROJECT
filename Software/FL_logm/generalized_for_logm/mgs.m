function [Q, R] = mgs(A,E)
%
%  function [Q, R] = mgs(A,E);
%
% modified Gram-Schmidt orthogonalization of the columns of
% A. The columns of A are assumed to be linearly independent.
%
% [Q, R] = mgs(A) returns a matrix Q with orthonormal columns
% and an invertible upper triangular matrix R so that A = Q*R.
% If only orthogonalization is wanted, the accumulation of R can be
% avoided by ommiting the output parameter:
% Q = mgs(A)
%
% [Q, R] = mgs(A,E) performes the orthonomalization in the E scalar
% product, i.e., E needs to be symmetric positive definite.
% 

  [m,n]=size(A);
  
  if (nargin==2)&&(any(any(E'-E)))
    error('MGS:input matrix E needs to be selfadjoint');
  end
  
  
  for k=1:n
    if nargout==2
      if nargin==2
        R(k,k) = Enorm(E,A(:,k));
      else
        R(k,k) = norm(A(:,k));
      end
      Q(:,k) = A(:,k)/R(k,k);
      for j=k+1:n
        if nargin==2
          R(k,j) = Q(:,k)'*(E*A(:,j));
        else
          R(k,j) = Q(:,k)'*A(:,j);
        end
        A(:,j) = A(:,j)-Q(:,k)*R(k,j);
      end
    else
      if nargin==2
        R = Enorm(E,A(:,k));
      else
        R = norm(A(:,k));
      end
      Q(:,k) = A(:,k)/R;
      for j=k+1:n
        if nargin==2
          R = Q(:,k)'*(E*A(:,j));
        else
          R = Q(:,k)'*A(:,j);
        end
        A(:,j) = A(:,j)-Q(:,k)*R;
      end
    end
  end
  
  function nrm=Enorm(E,x)
    nrm=sqrt(x'*(E*x));