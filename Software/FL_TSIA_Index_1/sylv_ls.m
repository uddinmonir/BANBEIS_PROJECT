function [ X, eigH ] = sylv_ls( A,H,M, Trans )
% function [ X, eigH ] = sylv_ls( A,H,M,GJ1,GJ2,J3,J4, Trans )
% function [ X, eigH ] = sylv_ls( A,H,M, Trans )
%
% Solve the Sylvester Equation 
% 
%      A*X + X*H + M = 0                                                (1)
% 
% with A a large and sparse and H small and dense matrices. 
% The algorithm is described in [1].
%
% Inputs: 
%  A        large and sparse input matrix of (1), dimension n0 x n0
%  H        small and dense input matrix of (1), dimension n1 x n1
%  M        right hand side of (1), dimension n0 x n1
%  Trans    Optional flag. If Trans=='T' 
%               A'*X + X*H' + M = 0                                     (2)
%           is solved instead of (1). 
%      
%
% Outputs: 
%  X        Solution of (1) or (2), truncated to a real matrix.
%  eigH     Eigenvalues of H
% 
%
% [1] Sparse-Dense Sylvester Equations in H2-Model Order Reduction;
%     Benner, Peter; Köhler, Martin; Saak, Jens;
%     MPI Magdeburg Preprints 2011. 
%
% Copyright 2011-2012, Martin Köhler
% MPI Magdeburg

if ( nargin < 4)
% if ( nargin < 8)
    Trans = 'N';
end

n0 = size(A,1); 
n1 = size(H,1);
if ( size(M,1) ~= n0 && size(M,2)~=n1) 
    error('The dimension of the matrix M does not fit.'); 
end

[U,S]=schur(H,'complex');
Mtilde = M*U;
Xtilde = sparse(n0,n1);
% I=speye(size(GJ1));
I=speye(size(A));

if (Trans=='N')
    for j = 1:n1
        rhs=- Mtilde(:,j);
        for i=1:j-1
            rhs = rhs - S(i,j)*Xtilde(:,i);
        end
         Xtilde(:,j)=(A+S(j,j)*I)\rhs;
  %         Xtilde(:,j)=[(GJ1+S(j,j)*I) GJ2;J3 J4]\rhs;
    end
else
    for jj=1:n1
        j=n1-jj+1;
        rhs=-Mtilde(:,j);
         for i=j+1:n1
            rhs = rhs - conj(S(j,i))*Xtilde(:,i);
         end
        Xtilde(:,j)=(A+S(j,j)*I)' \rhs;
    %      Xtilde(:,j)=[(GJ1+S(j,j)*I) GJ2;J3 J4]'\rhs;
    end
end
            
X=Xtilde*U';
if isreal(A)
X=real(X);
end
eigH=diag(S);
end

