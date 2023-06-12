function [Ar,Br,Cr,V,W] = IRKA_index1_for_logm(E1,J1,J2,J3,J4,B1,B2,C1,C2,r,MaxIter,Tol)
% function IRKA_index1_symmetric generate ROM
%  of the system:
%  E x'(t) = A x(t) + B u(t)   -----------(Ia)
%     y(t) = C x(t) + Da u(t)  -----------(Ib)
% where 
%
%   A = [ J1 J2 ]  E = [ M1 0 ]  B = [ B1 ]  and C = [ C1 C2 ]
%       [ J3 J4 ],     [ 0  0 ],     [ B2 ], 
%
%  A and E are symmetric and B=C^T
%
% ROM: (standard state-space model)
%
%  Er xr'(t) = Ar xr(t) + Br u(t)  ------(IIa)
%  yr(t)  = Cr x(t) + Dar u(t)  -------(IIb)
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
[n,k] = size(J2);
m = size(B1,2);
p = size(C1,1);
b = randn(m,r);
c = randn(p,r);
B=(B1-J2*(J4\B2));
C=(C1-C2*(J4\J3));
S = 100*rand(r,1);
bz=sparse(k,1);
cz=sparse(k,1);
%% Start iteration
for iter = 1:MaxIter
    S_old = S;
%% Compute projection subspaces
    V = zeros(n,r);
    W = zeros(n,r);
    j = 1;
    BT=B*b;
    CT=C'*c;
    while(j < r+1)        
      %  x = [S(j)*M1-J1 J2;J3 J4]\[B1*b(:,j);B2*b(:,j)];
        x = [S(j)*E1-J1 J2;J3 -J4]\[BT(:,j);bz];
        x = x(1:n,:);
        y = [S(j)*E1'-J1' J3';J2' -J4']\[CT(:,j);cz];
        y = y(1:n,:);
        if(abs(imag(S(j))) > 0)
            V(:,j) = real(x);
            W(:,j) = real(y);
            V(:,j+1) = imag(x);
            W(:,j+1) = imag(y);
            j = j + 2;
        else
            V(:,j) = real(x);
            W(:,j) = real(y);
            j = j + 1;
        end
    end    
    [V,~] = qr(V,0);
    [W,~] = qr(W,0);
%% Compute ROM
    Er = (W'*E1*V);
    J1til = W'*J1*V; 
    J2til = W'*J2; 
    J3til = J3*V;
    B1til = W'*B1; 
    C1til = C1*V;
    Ar = Er\(J1til-J2til*(J4\J3til));
    Br = Er\(B1til-J2til*(J4\B2));
    Cr = C1til-C2*(J4\J3til);
   %Dar = Da-C2*(J4\B2);
%% Update interpolation points/tangential directions
    [T,S] = eig(Ar);
    S = - diag(S);
    b = (T\(Er\Br)).';
    c = Cr*T;
%% Check for convergence
   % err = norm(sort(S)-sort(S_old))/norm(S_old);
%% comment out for non-verbose mode %%
  %  fprintf('IRKA step %d, conv. crit. = %e \n', iter, err)
%    if(err < Tol)        
  %      break
 %   end
end
%if (iter == MaxIter && err > Tol),
%  fprintf('IRKA: No convergence in %d iterations.\n', MaxIter)
%end