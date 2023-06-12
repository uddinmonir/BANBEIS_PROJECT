function [Ar,Br,Cr] = index1_TSIA(E1,J1,J2,J3,J4,B1,B2,C1,C2,Ar,Br,Cr,maxiter,tol)

%% Initialization
%------------------
[n,k] = size(J2);
m = size(B1,2);
p = size(C1,1);
E=E1;
A=J1-J2*(J4\J3);
B=B1-J2*(J4\B2);
C=C1-C2*(J4\J3);
bz=sparse(k,1);
cz=sparse(k,1);
%---------------
[T,S] = eig(Ar);
I=eye(size(A,1));
Ir=eye(size(Ar,1));
P = (E1\B);
Q = -C';
%% Start iteration
for iter = 1:maxiter
    S_old = S;
    % Sr_omega= (1i/2*pi)*logm((Ar+1i*w1*Ir)/(Ar-1i*w2*Ir));
    %  P = E1B*Br';
    %Q = -C'*Cr;
   % P = (S_omega*B*Br'+B*Br'*Ar_omg');
   % Q = -(S_omega'*C'*Cr-C'*Cr*Ar_omg);
 %  [ X, ~] = sylv_ls_index1(J1,J2,J3,J4,Ar',P*Br');
  % [ Y, ~] = sylv_ls_index1(J1',J3',J2',J4',Ar,Q*Cr);
  % [X, ~] = sylv_ls(A,Ar',P*Br');
  % [Y, ~] = sylv_ls(A',Ar,Q*Cr);
  % [X, ~] = sylv_ls_sg(A,E,Ar',P*Br');
  % [Y ,~] = sylv_ls_sg(A',E',Ar,Q*Cr);
  %[X ,~] = sylv_ls_sg_index1(E1,J1,J2,J3,J4,Ar',P*Br');
  %[Y ,~] = sylv_ls_sg_index1(E1',J1',J3',J2',J4',Ar,Q*Cr);
  [ X ,~] = index1_sylv(E1,J1,J2,J3,J4,Ar',P*Br');
   [ Y ,~] = index1_sylv(E1',J1',J3',J2',J4',Ar,Q*Cr);
% X = lyap(A,Ar',P*Br');
% Y = lyap(A',Ar,Q*Cr);
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
