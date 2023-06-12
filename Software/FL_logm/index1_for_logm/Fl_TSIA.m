function [Er,Ar,Br,Cr] = index1_Fl_TSIA(E,A,B,C,Er,Ar,Br,Cr,maxiter,tol,w)

%% Initialization
n = size(A,1);
m = size(B,2);
p = size(C,1);
%S = 100*rand(r,1);
% b = randn(m,r);
% c = randn(p,r);
% So = S; bo = b;co = c;
%     
% [Er, Ar, Br, Cr, S, V] = IRKA(E,A,B,C,r,maxiter,tol);
 S= diag(eig(Er\Ar));
%S = - diag(S);
%  b = (T\(Er\Br)).';
%  c = Cr*T;
%  w=w2;
%disp('Finished Conventional IRKA')
%Aexm=expm(A*Time);
I=eye(size(A,1));
Ir=eye(size(Ar,1));
%A_omg= real((1i/2*pi)*logm((A+1i*wom1*eye(size(A,1)))\(A+1i*wom2*eye(size(A,1)))));
A_omg= (1i/2*pi)*logm((A+1i*w*I)/(A-1i*w*I));
%AexmT=expm(A'*Time);
%% Start iteration
for iter = 1:maxiter
    S_old = S;
    %% Compute projection subspaces
    %V = zeros(n,r);
   % W = zeros(n,r);
    %j = 1;
    Ar_omg= (1i/2*pi)*logm((Ar+1i*w*Ir)/(Ar-1i*w*Ir));
   %B_omg= F_omg*(Br*Br')+(Br*Br')*F_omg';
    %X = B*b - Aexm*B*b*expm(-diag(S)*Time);
    P = (A_omg*B*Br'+B*Br'*Ar_omg');
    Q = -(A_omg'*C'*Cr-C'*Cr*Ar_omg);
    %
    X=lyap(A,Ar',P);
    Y= lyap(A',Ar,Q);
    [V,~] = qr(X,0);
    [W,~] = qr(Y,0);
    % 
    W=W/(V'*W);
%     rank(V);
%     rank(W);
%    pause(.1);
    %% Compute ROM
    
    Er = (W'*E*V);
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
