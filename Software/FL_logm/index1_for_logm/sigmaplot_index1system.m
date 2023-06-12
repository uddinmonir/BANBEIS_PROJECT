function [s,Ho,Hr,abserr,relaterr] =sigmaplot_index1system(E1,J1,J2,J3,J4,B1,B2,C1,C2,Er,Ar,Br,Cr,Dar,low,up,points)
%
%  function sigmaplot_2ndTo1st draw sigmaplot of original and reduced model along
%  with their absolulate and relative deviations within the range [low,up].
% Input:
%  
% 
%  Output:
%  -> s is a random vector on logspace in [low,up]
%  -> Ho, Hr, are vector of maximum Hankel singular values of original and
%         reduced models, i.e., Ho=||G(jw)||_inf and Hr= Hr=||Gr(jw)||_inf
%  -> abserr is a vector of maximum singularvalues of the difference of
%      original and reduced model transfer function, i.e
%      ||G(jw)-Gr(jw)||_inf
% -> relaterr = ||G(jw)-Gr(jw)||_inf/||G(jw)||_inf 
%
  [l1,l2] = size(J2);
     E = [E1 sparse(l1,l2); sparse(l2,l1) sparse(l2,l2)];
     A = [J1 J2;J3 J4];
     B = [B1;B2];
     C = [C1 C2];
    s=logspace(low,up, points);
      for k=1:points
          G1 = C(1,:)*((1j*s(k)*E-A)\B(:,1));
          G2 = Cr(1,:)*((1j*s(k)*Er-Ar)\Br(:,1))+Dar(1,1);
          Ho(k) = max(svds(G1));
          Hr(k) = max(svds(G2));
          abserr(k) = max(svds(G1-G2));
      end
          relaterr=abserr./Ho;
%  sigmaplot for full and ROM    
       figure(300);
       loglog(s,Ho,'r')
        hold on; 
       loglog(s,Hr,'b-.')
        xlabel('\omega')
        ylabel('\sigma_{max}(G(j\omega))')
        title('Transfer function of original and reduced order system ')
        legend('full model','ROM')
        hold off
%  absolute deviation
        figure(400);
        loglog(s,abserr,'b')
         xlabel('\omega')
         ylabel('\sigma_{max}(G(j\omega)-G_r(j\omega))')
         title('absolute model reduction error')
         hold off
%   relative error;
       figure(500);
       loglog(s,abserr./Ho,'b')
       xlabel('\omega')
       ylabel('\sigma_{max}(G(j\omega)-G_r(j\omega))/ \sigma_{max}(G(j\omega))')
       title('relative model reduction error')
   
end