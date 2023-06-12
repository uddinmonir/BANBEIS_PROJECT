function [Ar,V,W] = IRKA11(E,A,B,C,r,maxiter)
%% Initialization
n = size(A,1);
m = size(B,2);
p = size(C,1);
S = 100*rand(r,1);
b = randn(m,r);
c = randn(p,r);
So = S; bo = b;co = c;
%% Start iteration
for iter = 1:maxiter
    S_old = S;
    %% Compute projection subspaces
    V = zeros(n,r);
    W = zeros(n,r);
    j = 1;
    while(j < r+1)        
        x = (S(j)*E-A)\(B*b(:,j));
        y = (S(j)*E'-A')\(C'*c(:,j));
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
    %W=W/(V'*W);
    %% Compute ROM
    Er = (V'*E*V);
    Ar = Er\(V'*A*V);
    Br = Er\(V'*B);
    Cr = C*V;    
    %Er = eye(r);
   
    %% Update interpolation points/tangential directions
    [T,S] = eig(Ar);
    S = - diag(S);
    b = (T\(Br)).';
   % b = (T\(Er\Br)).';
    c = Cr*T;
    
    %% Check for convergence
   % err = norm(sort(S)-sort(S_old))/norm(S_old);
    %% comment out for non-verbose mode %%
  %  fprintf('IRKA step %d, conv. crit. = %e \n', iter, err)
  %  err_hist(iter) = err;
  %  if(err < tol)        
 %       break
  %  end
end
%if (iter == maxiter && err > tol),
 % fprintf('IRKA: No convergence in %d iterations.\n', maxiter)
%end