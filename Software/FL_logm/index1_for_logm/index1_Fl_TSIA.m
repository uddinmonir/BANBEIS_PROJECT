function [Ar,Br,Cr] = index1_Fl_TSIA(E1,J1,J2,J3,J4,B1,B2,C1,C2,Ar,Br,Cr,maxiter,tol,w1,w2)

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
%S = 100*rand(r,1);
%w=w2;
%b = randn(m,r);
%c = randn(p,r);
%So = S; bo = b;co = c;
%[Er, Ar, Br, Cr, S] = IRKA_index1(E1,J1,J2,J3,J4,B1,B2,C1,C2,r,maxiter,tol);
[T,S] = eig(Ar);
%S = - diag(S);
%  b = (T\(Er\Br)).';
tic;[Ar1,Br1,Cr1,V1,W1] = IRKA_index1_for_logm(E1,J1,J2,J3,J4,B1,B2,C1,C2,200,maxiter,tol);toc
%[Ar1,V1,W1] = IRKA11(E,A,B,C,150,maxiter);
Ir1=eye(size(Ar1,1));
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
%A_omg= real((1i/2*pi)*logm((A+1i*wom1*eye(size(A,1)))\(A+1i*wom2*eye(size(A,1)))));
%A_omg= (1i/2*pi)*logm((A+1i*w1*I)/(A-1i*w2*I));
%AexmT=expm(A'*Time);
%% Start iteration
for iter = 1:maxiter
    S_old = S;
    %% Compute projection subspaces
   Sr_omega= real((1i/pi)*logm((Ar+1i*w1*Ir)\(Ar+1i*w2*Ir)));
   %B_omg= F_omg*(Br*Br')+(Br*Br')*F_omg';
   % P = E*S_omega*(E\B)*Br'+B*Br'*Sr_omega';
   % Q = -(S_omega'*C'*Cr+C'*Cr*Sr_omega);
    %X = B*b - Aexm*B*b*expm(-diag(S)*Time);
    P = (S_omega*B*Br'+B*Br'*Sr_omega');
    Q = -(S_omega'*C'*Cr+C'*Cr*Sr_omega);
    %
     [ X ,~] = index1_sylv(E1,J1,J2,J3,J4,Ar',P);
     [ Y ,~] = index1_sylv(E1',J1',J3',J2',J4',Ar,Q);
   %[X ,~] = sylv_ls_sg_index1(E1,J1,J2,J3,J4,Ar',P);
   %[Y ,~] = sylv_ls_sg_index1(E1',J1',J3',J2',J4',Ar,Q);
   % [X ,~] = sylv_ls_sg(A,E,Ar',P);
   % [Y ,~] = sylv_ls_sg(A',E',Ar,Q);
   % X=lyap(A,Ar',P);
   % Y= lyap(A',Ar,Q);
    [V,~] = qr(X,0);
    [W,~] = qr(Y,0);
    % 
    W=W/(V'*W);
%     rank(V);
%     rank(W);
%    pause(.1);
    %% Compute ROM
    
    Er = (W'*E1*V);
    Ar = Er\(W'*A*V);
    Br = Er\(W'*B);
    Cr = C*V;    
    %% Update interpolation points/tangential directions
    S = diag(eig(Ar));
    %S = diag(S);
    %b = (T\(Er\Br)).';
    %c = Cr*T;
    
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
