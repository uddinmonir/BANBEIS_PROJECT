function [ntAr,ntBr,ntCr] = indexI_TSIA(A,B,C,Ar,Br,Cr,GJ1,GJ2,JJ43,JB42,GB1,C1,C2,J3,J4,Tmaxiter,tol)


S= eig(Ar);
ntAr = Ar;
ntBr = Br;
ntCr = Cr;
%% Start iteration
for iter = 1:Tmaxiter
    S_old = S;
    
        P = B*ntBr';
        Q = -(C' *ntCr);
   %    Q = (C' *ntCr);
       [ X ,~] = sylv_ls(A,ntAr',P);
       [ Y ,~] = sylv_ls(A',ntAr,Q);
%      [ X ,~] = sylv_ls(A,ntAr',P,GJ1,GJ2,J3,J4);
%      [ Y ,~] = sylv_ls(A',ntAr,Q,GJ1',J3',GJ2',J4');
%    X = sylvester(full(A),full(ntAr'),full(P));
%    Y = sylvester(full(A'),full(ntAr),full(Q));
 
    [V,~] = qr(X,0);
    [W,~] = qr(Y,0);
    % 
    W=W/(V' *W);
    %% Compute ROM
    ntAr = (W' *GJ1*V) - ((W' *GJ2)*(JJ43*V));
    ntBr =  (W' *GB1) - ((W' *GJ2)*JB42);
    ntCr = (C1*V) - (C2*(JJ43*V));

%% Update interpolation points/tangential directions
    S = eig(ntAr);
    %% Check for convergence
    err = norm(sort(S)-sort(S_old))/norm(S_old);
    %% comment out for non-verbose mode %%
    fprintf('TSIA step %d, conv. crit. = %e \n', iter, err)
   % err_hist(iter) = err;
    if(err < tol)        
        break
    end
end

if (iter == Tmaxiter && err > tol),
  fprintf('TSIA: No convergence in %d iterations.\n', Tmaxiter)
end
