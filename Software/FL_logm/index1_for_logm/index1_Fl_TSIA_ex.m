function [Ar,Br,Cr] = index1_Fl_TSIA_ex(E1,J1,J2,J3,J4,B1,B2,C1,C2,Ar,Br,Cr,maxiter,tol,w1,w2)

%% Initialization
[n,k] = size(J2);
m = size(B1,2);
p = size(C1,1);
% b = randn(m,r);
% c = randn(p,r);
E=E1;
A=J1-J2*(J4\J3);
B=B1-J2*(J4\B2);
C=C1-C2*(J4\J3);
bz=sparse(k,1);
cz=sparse(k,1);
%Ar=Er\Ar;
[T,S] = eig(Ar);
I=eye(size(A,1));
Ir=eye(size(Ar,1));
%A_omg= (1i/2*pi)*logm((A+1i*w1*I)/(A-1i*w2*I));
S_omega= (A+1i*w1*E)\(A+1i*w2*E);
%% Patrick
S_omega = real((1i/pi)*logm(full(S_omega)));
%% Start iteration
for iter = 1:maxiter
    S_old = S;
    Sr_omega= (1i/2*pi)*logm((Ar+1i*w1*Ir)/(Ar-1i*w2*Ir));
     P = E*S_omega*(E\B)*Br'+B*Br'*Sr_omega';
    Q = -(S_omega'*C'*Cr+C'*Cr*Sr_omega);
   % P = (S_omega*B*Br'+B*Br'*Ar_omg');
   % Q = -(S_omega'*C'*Cr-C'*Cr*Ar_omg);
   [ X, ~] = sylv_ls_index1(J1,J2,J3,J4,Ar',P);
   [ Y, ~] = sylv_ls_index1(J1',J3',J2',J4',Ar,Q);
%[ X, ~] = sylv_ls(A,Ar',P);
%[ Y, ~] = sylv_ls( A',Ar,Q);
    [V,~] = qr(X,0);
    [W,~] = qr(Y,0);
    % 
    W=W/(V'*W);
    %% Compute ROM
    Er = (W'*E*V);
    Ar = Er\(W'*A*V);
    Br = Er\(W'*B);
    Cr = C*V;  
    %% Update interpolation points/tangential directions
    S = diag(eig(Ar));
    
    %% Check for convergence
    err = norm(sort(S)-sort(S_old))/norm(S_old);
    %% comment out for non-verbose mode %%
    fprintf('IRKA step %d, conv. crit. = %e \n', iter, err)
    err_hist(iter) = err;
    if(err < tol)        
        break
    end
end

if (iter == maxiter && err > tol),
  fprintf('IRKA: No convergence in %d iterations.\n', maxiter)
end
