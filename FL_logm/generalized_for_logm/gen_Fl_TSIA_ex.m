function [Ar,Br,Cr] = gen_Fl_TSIA_ex(E,A,B,C,Ar,Br,Cr,V,W,maxiter,tol,w1,w2)

%% Initialization
n = size(A,1);
m = size(B,2);
p = size(C,1);
S= eig(Ar);
I=eye(size(A,1));
Ir=eye(size(Ar,1));
[Ar1,V1,W1] = IRKA_gen_for_logm(E,A,B,C,150,maxiter);
Ir1=eye(size(Ar1,1));
%Ar = IRKA(E,A,B,C,r,maxiter,tol)
%[V]=rksm(E,A,B,4000,100,10^(-15),10^(-15));
%S_omega= (A+1i*w1*E)\(A+1i*w2*E);
%Atil=W'*A*V;  Etil=W'*E*V;
%S_omega= (Atil+1i*w1*Etil)\(Atil+1i*w2*Etil);
S_omega= real((1i/pi)*V1*logm((Ar1+1i*w1*Ir1)\(Ar1+1i*w2*Ir1))*V1');
%% Start iteration
for iter = 1:maxiter
    S_old = S;
    %% Patrick
 Sr_omega= real((1i/pi)*logm((Ar+1i*w1*Ir)\(Ar+1i*w2*Ir)));
 
 P = E*S_omega*(E\B)*Br'+B*Br'*Sr_omega';
 %P = S_omega*B*Br'+B*Br'*Sr_omega';
 Q = -(S_omega'*C'*Cr+C'*Cr*Sr_omega);
    % danse-sparse matrix equation
    %[ X, ~] = sylv_ls(A,Ar',P);
    [ X ,~] = sylv_ls_sg(A,E,Ar',P);
    [ Y ,~] = sylv_ls_sg(A',E',Ar,Q);
    %Y= lyap(A',Ar,-Q);
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
    S = eig(Ar);
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
