function [Art,Brt,Crt] = indexI_tl_TSIA(A,B,C,Ar,Br,Cr,GJ1,GJ2,JJ43,JB42,GB1,C1,C2,J3,J4,Tmaxiter,tol,T1,T2,Ctol,nZtol)

S=eig(Ar);
Art = Ar;
Brt = Br;
Crt = Cr;
Arexm = customExpm((A*T2),Ctol,nZtol);
ArexmT = customExpm((A' *T2),Ctol,nZtol);
%% Start iteration
for iter = 1:Tmaxiter
    S_old = S;
    Artexm = customExpm((Art*T2),Ctol,nZtol);
    ArtexmT = customExpm((Art' *T2),Ctol,nZtol);
   
    P = (B*Brt') - ((Arexm*B)*(Brt' *ArtexmT));
    Q = ((ArexmT*C') *(Crt*Artexm)) - (C' *Crt);
%   Q = (C' *Crt) - ((AexmT*C') *(Crt*Artexm)) ;
    
    % dense-sparse matrix equation
     [ X ,~] = sylv_ls(A,Art',P);
     [ Y ,~] = sylv_ls(A',Art,Q);
%         [ X ,~] = sylv_ls(A,Art',P,GJ1,GJ2,J3,J4);
%        [ Y ,~] = sylv_ls(A',Art,Q,GJ1',J3',GJ2',J4');
%    X = sylvester(full(A),full(Art'),full(P));
%    Y = sylvester(full(A'),full(Art),full(Q));
    

    [V,~] = qr(X,0);
    [W,~] = qr(Y,0);
    % 
    W=W/(V' *W);
    %% Compute ROM
    Art = (W' *GJ1*V) - ((W' *GJ2)*(JJ43*V));
    Brt = (W' *GB1) - ((W' *GJ2)*JB42);
    Crt = (C1*V) - (C2*(JJ43*V));    

    %% Update interpolation points/tangential directions
    %S=eig(Art);
    %% Check for convergence
    %err = norm(sort(S)-sort(S_old))/norm(S_old);
    %% comment out for non-verbose mode %%
%     fprintf('TL-TSIA step %d, conv. crit. = %e \n', iter, err)
      fprintf('TL-TSIA step %d \n', iter)
   % err_hist(iter) = err;
%     if(err < tol)        
%         break
%     end
end

% if (iter == Tmaxiter && err > tol),
%   fprintf('TL-TSIA: No convergence in %d iterations.\n', Tmaxiter)
% end
