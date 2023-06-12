function [Er,Ar,Br,Cr] = fl_irka_generalized(E,A,B,C,r,S,Er,Ar,b,c,maxiter,tol,wom1,wom2)

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
% [T,S] = eig(Ar,Er);
%S = - diag(S);
% b = (T\(Er\Br)).';
% c = Cr*T;
%disp('Finished Conventional IRKA')
%Aexm=expm(A*Time);
A_omg= real((1i/pi)*log((A+1i*wom1*E)\(A+1i*wom2*E)));
%AexmT=expm(A'*Time);
%% Start iteration
for iter = 1:maxiter
    S_old = S;
    %% Compute projection subspaces
    V = zeros(n,r);
    W = zeros(n,r);
    j = 1;
    % X = B*b;
    % Y = C'*c;
      Ar_omg= real((1i/pi)*log((Ar+1i*wom1*Er)\(Ar+1i*wom2*Er)));
   %B_omg= F_omg*(Br*Br')+(Br*Br')*F_omg';
   %Y=lyap(Ar,B_omg);  
    %X = B*b - Aexm*B*b*expm(-diag(S)*Time);
    X = A_omg*B*b+B*b*Ar_omg';
    %Y = C'*c - AexmT*C'*c*expm(-diag(S)*Time); 
   Y = A_omg'*C'*c+C'*c*Ar_omg;
    while(j < r+1)        
        x = (S(j)*E-A)\(X(:,j)); 
        y = (S(j)*E'-A')\(Y(:,j));
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
%     rank(V);
%     rank(W);
%    pause(.1);
    %% Compute ROM
    Er = (W'*E*V);
    Ar = Er\(W'*A*V);
    Br = Er\(W'*B);
    Cr = C*V;    
    Er = eye(r);
    %% Update interpolation points/tangential directions
    [T,S] = eig(Ar,Er);
    S = - diag(S);
    b = (T\(Er\Br)).';
    c = Cr*T;
    
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
