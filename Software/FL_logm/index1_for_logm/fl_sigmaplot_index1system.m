function [s,tf] =fl_sigmaplot_index1system(E1,J1,J2,J3,J4,B1,B2,C1,C2,Er,Ar,Br,Cr,Er_fl,Ar_fl,Br_fl,Cr_fl,Dar,low,up,points)
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
    s=logspace(low,up,points);
      for k=1:points
          G1 = C*((1j*s(k)*E-A)\B);
          G2 = Cr*((1j*s(k)*Er-Ar)\Br)+Dar;
          G3=Cr_fl*((1j*s(k)*Er_fl-Ar_fl)\Br_fl)+Dar;

          tf(k)=max(svds(G1));
          tf_rom(k)=max(svds(G2));
          tf_rom_fl(k) = max(svds(G3));
          abs_err(k)=max(svds(G1-G2));
          abs_err_fl(k)=max(svds(G1-G3));
          
      end
          rel_err=abs_err./tf;
          rel_err_fl=abs_err_fl./tf;
%  Transfer Function Plot for full,frequency limited and frequency unrestricted reduced model    
        figure(10);
        loglog(s,tf,'k')
        hold on;
        loglog(s,tf_rom,'b--')
        hold on;
        loglog(s,tf_rom_fl,'r-.')
        xlabel('\omega')
        ylabel('\sigma_{max}(G(j\omega))')
        title('Transfer function')
        legend('Full Model','Frequency Unrestricted', 'Frequency Restricted')
        hold off
%  Absolute Error
        figure(20);
       loglog(s,abs_err,'b--')
        hold on;
        loglog(s, abs_err_fl, 'r-.')
         xlabel('\omega')
         ylabel('\sigma_{max}(G(j\omega)-G_r(j\omega))')
         title('absolute error')
         legend('Frequency unrestricted', 'Frequency Restricted')
         hold off
%   relative error;
       figure(30);
       loglog(s, rel_err,'b--')
       hold on;
       loglog(s,rel_err_fl,'r-.')
       xlabel('\omega')
       ylabel('\sigma_{max}(G(j\omega)-G_r(j\omega))/ \sigma_{max}(G(j\omega))')
       title('Relative error')
       legend('Frequency Unrestricted', 'Frequency Restricted')
 hold off
end