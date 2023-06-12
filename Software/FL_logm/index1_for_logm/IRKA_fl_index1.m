function [Ar,Br,Cr,S] = IRKA_fl_index1(E1,J1,J2,J3,J4,B1,B2,C1,C2,Ar, Br, Cr,r,maxiter,tol,w1,w2)

%% Initialization
[n,k] = size(J2);
m = size(B1,2);
p = size(C1,1);
% b = randn(m,r);
% c = randn(p,r);
A=J1-J2*(J4\J3);
B=B1-J2*(J4\B2);
C=C1-C2*(J4\J3);
bz=sparse(k,1);
cz=sparse(k,1);
[T,S] = eig(Ar);
%S = - diag(S);
b = (T\(Br)).';
c = Cr*T;
S = - diag(S);
%disp('Finished Conventional IRKA')
A_omg= real(1i/pi)*logm(full((A+1i*w1*E1)\(A+1i*w2*E1)));
%A_omg= (1i/2*pi*logm(-A-1i*w*eye(size(A))));
%Exp_A=expm(A*Time)*B;
%TExp_A=expm(A'*Time)*C';
%% Start iteration
for iter = 1:maxiter
    S_old = S;
    %% Compute projection subspaces
    V = zeros(n,r);
    W = zeros(n,r);
    j = 1;
    % X = B*b;
    Ar_omg= real((1i/pi)*logm((Ar+1i*w1*eye(size(Ar)))\(Ar+1i*w2*eye(size(Ar)))));
    %Ar_omg= real(1i/pi*logm(-(Er\Ar)-1i*w*eye(size(Ar))));
    % Y = C'*c;
    X = A_omg*B*b+B*b*Ar_omg';
    %Y = C'*c - AexmT*C'*c*expm(-diag(S)*Time); 
   Y = (A_omg'*C'*c+C'*c*Ar_omg);
    % X = B*b - Exp_A*b*expm(-diag(S)*Time);
    % Y = C'*c - TExp_A*c*expm(-diag(S)*Time); 
    while(j < r+1)        
        x = [S(j)*E1-J1 -J2;-J3 -J4]\[X(:,j);bz];
        x = x(1:n,:);
        y = [S(j)*E1'-J1' -J3';-J2' -J4']\[Y(:,j);cz];
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
    [V,~] = qr(V,0);
    [W,~] = qr(W,0);
   pause(.1)
    %% Compute ROM
    Er = (W'*E1*V);
    J1til = W'*J1*V; 
    J2til = W'*J2; 
    J3til = J3*V;
    B1til = W'*B1; 
    C1til = C1*V;
    Ar = Er\(J1til-J2til*(J4\J3til));
    Br = Er\(B1til-J2til*(J4\B2));
    Cr = C1til-C2*(J4\J3til);
    %% Update interpolation points/tangential directions
    [T,S] = eig(Ar);
    S = - diag(S);
    b = (T\(Br)).';
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
