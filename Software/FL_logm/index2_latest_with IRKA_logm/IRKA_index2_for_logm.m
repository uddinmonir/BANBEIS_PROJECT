function  [Ar,Br,V] = IRKA_index2_for_logm(M,A1,A2,B1,r,MaxIter,Tol)
% function IRKA_index1_symmetric generate ROM
%  of the system:
%  E x'(t) = A x(t) + B u(t)   -----------(Ia)
%     y(t) = C x(t) + Da u(t)  -----------(Ib)
% where 
%
%   A = [ J1 J2 ]  E = [ M1 0 ]  B = [ B1 ]  and C = [ C1 C2 ]
%       [ J3 0 ],     [ 0  0 ],     [ B2 ], 
%
%  A and E are symmetric and B=C^T
%
% ROM: (standard state-space model)
%
%  Er xr'(t) = Ar xr(t) + Br u(t)  ------(IIa)
%     yr(t) = Cr x(t) + Dar u(t)  -------(IIb)
%
% INPUT :
%     M1, J1, J2, J3, J4, B1, B2, C1, C2, Da matrices as in (I) 
%
% OUTPUT:
%     Er, Ar, Br, Cr, and Dar as in (II)
%
%  Mohammad Monir Uddin, January 2014, MPI Magdeburg
%
%% Initialization
[n1,n2] = size(A2);
m = size(B1,2);
S = 100*rand(r,1);
%p = size(C1,1);
b = randn(m,r);
%c = randn(p,r);
A4=sparse(n2,n2);
B2=sparse(n2,m);
%C2=sparse(p,n2);
%% Start iteration
for iter = 1:MaxIter
    S_old = S;
%% Compute projection subspaces
    V = zeros(n1,r);
    %W = zeros(n1,r);
    j = 1;
    while(j < r+1)        
      %  x = [S(j)*M1-J1 J2;J3 J4]\[B1*b(:,j);B2*b(:,j)];
        x = [A1-S(j)*M A2;A2' A4]\[B1*b(:,j);B2*b(:,j)];
        x = x(1:n1,:);
%         y = [A1'-S(j)*M' A2;A2' A4]\[C1'*c(:,j);C2'*c(:,j)];
%         y = y(1:n1,:);
        if(abs(imag(S(j))) > 0)
            V(:,j) = real(x);
            %W(:,j) = real(y);
            V(:,j+1) = imag(x);
           % W(:,j+1) = imag(y);
            j = j + 2;
        else
            V(:,j) = real(x);
            %W(:,j) = real(y);
            j = j + 1;
        end
    end    
    [V,~] = qr((V),0);
   % [W,~] = qr((W),0);
%% Compute ROM
    Er = (V'*M*V);
    Ar = Er\(V'*A1*V); 
    Br = Er\(V'*B1); 
   % Cr = C1*V;
    %% Update interpolation points/tangential directions
    [T,S] = eig(Ar);
    S = - diag(S);
    b = (T\Br).';
   % c = Cr*T;
%% Check for convergence
%     err = norm(sort(S)-sort(S_old))/norm(S_old);
% %% comment out for non-verbose mode %%
%     fprintf('IRKA step %d, conv. crit. = %e \n', iter, err)
%     if(err < Tol)        
%         break
%     end
% end
% if (iter == MaxIter && err > Tol),
%   fprintf('IRKA: No convergence in %d iterations.\n', MaxIter)
end