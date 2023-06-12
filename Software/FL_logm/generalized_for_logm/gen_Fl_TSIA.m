function [Ar,Br,Cr] = gen_Fl_TSIA(E,A,B,C,Ar,Br,Cr,maxiter,tol,w1,w2)

%% Initialization
n = size(A,1);
m = size(B,2);
p = size(C,1);
S= eig(Ar);
I=eye(size(A,1));
Ir=eye(size(Ar,1));
%[V]=rksm(E,A,B,4000,100,10^(-15),10^(-15));
S_omega= (A+1i*w1*E)\(A+1i*w2*E);
%Atil=V'*A*V;  Etil=V'*E*V;
%S_omega= (Atil+1i*w1*Etil)\(Atil+1i*w2*Etil);
%S_omega= -logm_Custom(A+1i*w1*E)+logm_Custom(A+1i*w2*E);
%% Patrick
S_omega = real((1i/pi)*logm((full(S_omega))));
%S_omega = real((1i/pi)*logm_Custom((S_omega)));
%LM1=real((1i/pi)*logm_Custom(A+1i*w1*E));
%LM2=real((1i/pi)*logm_Custom(A+1i*w2*E));
%S_omega = real((1i/pi)*(logm_Custom(A+1i*w1*E)-logm_Custom(A+1i*w2*E)));
%S_omega=LM1-LM2;
%% Start iteration
for iter = 1:maxiter
    S_old = S;
    %% Patrick
 Sr_omega= real((1i/pi)*logm((Ar+1i*w1*Ir)\(Ar+1i*w2*Ir)));
 P = E*S_omega*(E\B)*Br'+B*Br'*Sr_omega';
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
