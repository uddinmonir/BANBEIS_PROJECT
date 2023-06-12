function [Ar,Br,Cr] = index2_Fl_TSIA(E1,A1,A2,B1,C1,Ar,Br,Cr,maxiter,tol,w1,w2)

%% Initialization
n1=size(A2,2);
A4=sparse(n1,n1);
I=speye(size(A1,1));
PI=I-A2*((A2'*(E1\A2))\(A2'/E1));
n = size(A1,1);
m = size(B1,2);
p = size(C1,1);
S= eig(Ar);
Ir=eye(size(Ar,1));
At=PI*A1;Et=E1;Bt=PI*B1;Ct=C1;
%[V]=rksm(E,A,B,4000,100,10^(-15),10^(-15));
S_omega= (At+1i*w1*E1)\(At+1i*w2*E1);
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
 %P = S_omega*(B*Br')+(B*Br')*Sr_omega';%
 P = E1*S_omega*(E1\Bt)*Br'+Bt*Br'*Sr_omega';
 Q = -(S_omega'*C1'*Cr+C1'*Cr*Sr_omega);
    % danse-sparse matrix equation
    [ X ,~] = index2_sylv(E1,A1,A2,A4,Ar',P);
    [ Y ,~] = index2_sylv(E1',A1',A2,A4,Ar,Q);
    %[ X ,~] = sylv_ls_sg(At,E1,Ar',P);
    %[ Y ,~] = sylv_ls_sg(At',E1',Ar,Q);
   % X = lyap(E\A,Ar',P);
   % Y = lyap((E\A)',Ar,Q);
    %Y= lyap(A',Ar,-Q);
    [V,~] = qr(X,0);
    [W,~] = qr(Y,0);
    % 
    W=W/(V'*W);
    %% Compute ROM
    
    Er = (W'*E1*V);
    Ar = Er\(W'*A1*V);
    Br = Er\(W'*B1);
    Cr = C1*V;    
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
