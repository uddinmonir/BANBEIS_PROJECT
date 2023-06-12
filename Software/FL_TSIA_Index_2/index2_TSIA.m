function [Ar,Br,Cr,V,W] = index2_TSIA(M,A1,A2,B1,C1,Ar,Br,Cr,MaxIter,Tol)

%% Initialization
%------------------
[~,n2] = size(A2);
%A3 = A2';
A4 = sparse(n2,n2);
%---------------
S = eig(Ar);
%P = B1*Br';
%Q = -C1'*Cr;
%% Start iteration
for iter = 1:MaxIter
    S_old=S;
   P=B1*Br';
   Q=-C1'*Cr;
    [ X ,~] = index2_sylv(M,A1,A2,A4,Ar',P);
    [ Y ,~] = index2_sylv(M',A1',A2,A4',Ar,Q);
    % X = lyap(A1,Ar',P*Br');
    % Y = lyap(A1',Ar,Q*Cr);
    [V,~] = qr(X,0);
    [W,~] = qr(Y,0);
    
    %% Compute ROM
    Er = (W'*M*V);
    Ar = Er\(W'*A1*V); 
    Br = Er\(W'*B1); 
    Cr = C1*V;
    %% Update interpolation points/tangential directions
    S = eig(Ar);
    
    %% Check for convergence
    err = norm(sort(S)-sort(S_old))/norm(S_old);
    %% comment out for non-verbose mode %%
    fprintf('TSIA step %d, conv. crit. = %e \n', iter, err)
    err_hist(iter) = err;
    if(err < Tol)
        break
    end
end

if (iter == MaxIter && err > Tol)
    fprintf('TSIA: No convergence in %d iterations.\n', MaxIter)
end
