function [ X ,eigH] = sylv_ls_sg(A,E,H,M,Trans )
% function [ X, eigH ] = sylv_ls_sg( A,E,H,M, Trans )
%
% Solve the semi-generalized Sylvester Equation 
% 
%      A*X + E*X*H + M = 0                                              (1)
% 
% with A and E large and sparse and H small and dense matrices. 
% The algorithm is described in [1].
%
% Inputs: 
%  A,E      large and sparse input matrix of (1), dimension n0 x n0
%  H        small and dense input matrix of (1), dimension n1 x n1
%  M        right hand side of (1), dimension n0 x n1
%  Trans    Optional flag. If Trans=='T' 
%               A'*X + E'*X*H' + M = 0                                  (2)
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

if ( nargin < 5)
    Trans = 'N';
end

n0 = size(A,1); 
n1 = size(H,1);
if ( size(M,1) ~= n0 && size(M,2)~=n1) 
    error('The dimension of the matrix M does not fit.'); 
end
if (size(A) ~= size(E)) 
    error('Dimension of A and E does not fit.');
end

q = n1;

[U,S]=schur(H,'complex');
Mtilde = M*U;
Xtilde = zeros(size(A,1),q);

if (Trans=='N')
    for j = 1:q
%         rhs=-Mtilde(:,j);
        rhs = zeros(size(Mtilde,1),1);
        for i=1:j-1
            rhs = rhs + S(i,j)*Xtilde(:,i);
        end
        rhs = -Mtilde(:,j) - E*rhs;
        Xtilde(:,j)=(A+S(j,j)*E)\rhs;
   
    end
else
    %S=S';
    for jj=1:q
        j=q-jj+1;
%         rhs=-Mtilde(:,j);
        rhs = zeros(size(Mtilde,1),1);
        for i=j+1:q
            rhs = rhs + conj(S(j,i))*Xtilde(:,i);
        end
        rhs = -Mtilde(:,j) - E'*rhs;
        Xtilde(:,j)=(A+S(j,j)*E)'\rhs;
    end
end
            
X=Xtilde*U';
X=real(X);
eigH=diag(S);
end

