function [Ar,Br,Cr] = index2_Fl_irka(E1,A1,A2,B1,C1,Ar,Br,Cr,r1,maxiter,tol,w1,w2)

%% Initialization
r=size(Ar,1);
[n1,n2] = size(A2);
m = size(B1,2);
S = 100*rand(r1,1);
p = size(C1,1);
b = randn(m,r);
c = randn(p,r);
A4=sparse(n2,n2);
B2=sparse(n2,1);
%C2=sparse(,r);
I=speye(size(A1,1));
PI=I-A2*((A2'*(E1\A2))\(A2'/E1));
n = size(A1,1);
m = size(B1,2);
p = size(C1,1);
S= eig(Ar);
Ir=eye(size(Ar,1));
At=PI*A1;Et=E1;Bt=PI*B1;Ct=C1;
[Ar1,Br1,V1] = IRKA_index2_for_logm(E1,A1,A2,B1,r1,maxiter,tol);
%[Ar1,V1] = IRKA_gen_for_logm(Et,At,Bt,150,maxiter);
Ir1=eye(size(Ar1,1));
%Ar = IRKA(E,A,B,C,r,maxiter,tol)
%[V]=rksm(E,A,B,4000,100,10^(-15),10^(-15));
%S_omega= (A+1i*w1*E)\(A+1i*w2*E);
%Atil=W'*A*V;  Etil=W'*E*V;
%S_omega= (Atil+1i*w1*Etil)\(Atil+1i*w2*Etil);
%S_omega= real((1i/pi)*V1*logm((Ar1+1i*w1*Ir1)\(Ar1+1i*w2*Ir1))*V1');
%Ir1=eye(size(Ar1,1));
%Ar = IRKA(E,A,B,C,r,maxiter,tol)
%[V]=rksm(E,A,B,4000,100,10^(-15),10^(-15));
S_omega= (Ar1+1i*w1*Ir1)\(Ar1+1i*w2*Ir1);
%Atil=W'*A*V;  Etil=W'*E*V;
%S_omega= (Atil+1i*w1*Etil)\(Atil+1i*w2*Etil);
%S_omega= real((1i/pi)*V1*logm((Ar1+1i*w1*Ir1)\(Ar1+1i*w2*Ir1))*V1');
%  c = Cr*T;
%  w=w2;
%disp('Finished Conventional IRKA')
%Aexm=expm(A*Time);
% I=eye(size(A,1));
 Ir=eye(size(Ar,1));
% S_omega= (A+1i*w1*E)\(A+1i*w2*E);
%% Patrick
S_omega = real((1i/pi)*V1*logm(S_omega)*V1');
%Atil=V'*A*V;  Etil=V'*E*V;
%S_omega= (Atil+1i*w1*Etil)\(Atil+1i*w2*Etil);
%S_omega= -logm_Custom(A+1i*w1*E)+logm_Custom(A+1i*w2*E);
%% Patrick
V = zeros(n,r);
W = zeros(n,r);
%S_omega = real((1i/pi)*logm((full(S_omega))));
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
 X = E1*S_omega*(E1\Bt)*Br'+Bt*Br'*Sr_omega';
 size(X)
 Y = -(S_omega'*C1'*Cr+C1'*Cr*Sr_omega);
  size(X)
 j = 1;
    % danse-sparse matrix equation
   % [ X ,~] = index2_sylv(E1,A1,A2,A4,Ar',P);
    %[ Y ,~] = index2_sylv(E1',A1',A2,A4,Ar,Q);
     while(j < r+1)        
      %  x = [S(j)*M1-J1 J2;J3 J4]\[B1*b(:,j);B2*b(:,j)];
        x = [A1-S(j)*E1 A2;A2' A4]\[X(:,j);B2];
        x = x(1:n1,:);
        y = [A1'-S(j)*E1' A2;A2' A4]\[Y(:,j);B2];
        y = y(1:n1,:);
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
    [T,S] = eig(Ar);
    S = - diag(S);
    b = (T\(Er\Br)).';
     c = Cr*T;
    %% Update interpolation points/tangential directions
    %S = eig(Ar);
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
