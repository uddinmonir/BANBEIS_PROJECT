function [Ar,Br,Cr] = gen_TSIA(E,A,B,C,Ar,Br,Cr,maxiter,tol)

%% Initialization
n = size(A,1);
m = size(B,2);
p = size(C,1);
S= eig(Ar);
I=eye(size(A,1));
Ir=eye(size(Ar,1));
%% Start iteration
for iter = 1:maxiter
    S_old = S;
    P = B*Br';
    Q = -C'*Cr;
    [ X ,~] = sylv_ls_sg(A,E,Ar',P);
    [ Y ,~] = sylv_ls_sg(A',E',Ar,Q);
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
