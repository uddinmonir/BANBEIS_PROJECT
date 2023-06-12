function [X,M,k] = sqrtm_DenBea(A,maxit)

n = length(A);
I = speye(n);
tol = 1e-30;
X = A;
M = A;
if nargin<2
maxit = 150;
end
for k = 1:maxit
    
   invM = I/M;
   X = X*(I + invM)/2;
   M = 0.5*(I + (M + invM)/2);
   Mres = norm(M - I,'fro');
   if Mres <= tol
       return; 
   end
fprintf('Iteranation No. %d \n', k)
end

