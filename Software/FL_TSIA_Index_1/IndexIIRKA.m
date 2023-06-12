function [Ar,Br,Cr] = IndexIIRKA(GJ1,GJ2,JJ43,GB1,JB42,J3,J4,B2,C1,C2,r,Imaxiter,tol,n,m,p)

%% Initialization

S = 100*rand(r,1);
b = randn(m,r);
c = randn(p,r);
I = speye(size(GJ1));   

%% Start iteration
for iter = 1:Imaxiter
    S_old = S;
    %% Compute projection subspaces
    V = zeros(n,r);
    W = zeros(n,r);
    j = 1; 
    while(j < r+1)        
        x = [((S(j)*I)-GJ1) GJ2;J3 J4]\[(GB1*b(:,j));(B2*b(:,j))];
        y = [((S(j)*I)-GJ1)' J3';GJ2' J4']\[(C1'*c(:,j));(C2'*c(:,j))];
        x = x(1:n,:);
        y = y(1:n,:);
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
    [Tr,~] = qr(V,0);
    [Tl,~] = qr(W,0);
  
    %% Compute ROM
    Ar = (Tl' *GJ1*Tr) - ((Tl' *GJ2)*(JJ43*Tr));        
    Br = (Tl' *GB1) - ((Tl' *GJ2)*JB42);             
    Cr = (C1*Tr) - (C2*(JJ43 *Tr));   

    %% Update interpolation points/tangential directions
    [T,S] = eig(Ar);
    S = - diag(S);
    b = (T\Br)';            
    c = Cr*T;
    
    %% Check for convergence
    err = norm(sort(S)-sort(S_old))/norm(S_old);
    %% comment out for non-verbose mode %%
    fprintf('IRKA step %d, conv. crit. = %e \n', iter, err)
   % err_hist(iter) = err;
    if(err < tol)        
        break
    end
end
if (iter == Imaxiter && err > tol),
  fprintf('IRKA: No convergence in %d iterations.\n', Imaxiter)
end