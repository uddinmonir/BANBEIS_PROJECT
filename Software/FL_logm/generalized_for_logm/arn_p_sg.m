function [H,V] = arn_p_sg(A,E,st, k,rv)
%
%  Arnoldi method w.r.t. (A,E) with index 1
%
%  Calling sequence:
%
%    [H,V] = arn_p(A,E,st,k)
%    [H,V] = arn_p(A,E,st,k,rv)
%
%  Input:
%
%    A,E       Input matrices for which the Arnoldi Algorithm is supposed to run.
%   st         number of states of descriptor system
%    k         number of Arnoldi steps (usually k << n);
%    rv        initial n-vector 
%              (optional - chosen by random, if omitted).
%    a         Alpha-shift  
%
%  Output:
%
%    H         matrix H ((k+1)-x-k matrix, upper Hessenberg);
%    V         matrix V (n-x-(k+1) matrix, orthogonal columns).
%
%  Method:
%
%    The Arnoldi method produces matrices V and H such that
%
%      V(:,1) in span{r},
%      V'*V = eye(k+1),
%      F*V(:,1:k) = V*H.
%
%  Remark:
%
%    This implementation does not check for (near-)breakdown!
%
   

% Input data not completely checked!

na = nargin;

n = size(A,1);                 % Get system order.


if k >= st-1, error('k must be smaller than the order of A!'); end
if na<5 || isempty(rv), rv = randn(st,1); end 

V = zeros(st,k+1);
H = zeros(k+1,k);

V(:,1) = (1.0/norm(rv))*rv;

beta = 0;

for j = 1:k
 
  if j > 1
    H(j,j-1) = beta;
    V(:,j) = (1.0/beta)*rv;
  end
  
  if (mod(j,5)==0)
    V(:,1:j) = mgs(V(:,1:j));
  end 

  w = E(1:st,1:st)\(A(1:st,1:st)*V(:,j));
  w = w - A(1:st,st+1:n)*(A(st+1:n,st+1:n)\(A(st+1:n,1:st)*V(:,j)));
  
  rv = w;
  
  for i = 1:j
    H(i,j) = V(:,i)'*w;
    rv = rv-H(i,j)*V(:,i);
  end

  beta = norm(rv);
  H(j+1,j) = beta;
 
end  

V(:,k+1) = (1.0/beta)*rv;