function [Ar,Br,Cr] = index2_Fl_TSIA(E1,A1,A2,B1,C1,Ar,Br,Cr,r1,maxiter,tol,w1,w2)

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
[Ar1,~,V1] = IRKA_index2_for_logm(E1,A1,A2,B1,r1,maxiter,tol);
%[Ar1,V1] = IRKA_gen_for_logm(Et,At,Bt,150,maxiter);
Ir1=eye(size(Ar1,1));
S_omega= (Ar1+1i*w1*Ir1)\(Ar1+1i*w2*Ir1);
 Ir=eye(size(Ar,1));
% S_omega= (A+1i*w1*E)\(A+1i*w2*E);
%% Patrick
S_omega = real((1i/pi)*V1*logm(S_omega)*V1');
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
   % X = lyap(Et\At,Ar',P);
   % Y = lyap((Et\At)',Ar,Q);
    %Y= lyap(A',Ar,-Q);
  % X = sylvester(Et\At,Ar',P);
   %Y = sylvester((Et\At)',Ar,Q);
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
