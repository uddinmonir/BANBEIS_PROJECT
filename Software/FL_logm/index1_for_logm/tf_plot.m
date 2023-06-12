function [space,tf] =tf_plot(E1,J1,J2,J3,J4,B1,B2,C1,C2,Er,Ar,Br,Cr,Er_fl,Ar_fl,Br_fl,Cr_fl,Dar,low_point,up_point,tot_point)
     
[l1,l2] = size(J2);
     E = [E1 sparse(l1,l2); sparse(l2,l1) sparse(l2,l2)];
     A = [J1 J2;J3 J4];
     B = [B1;B2];
     C = [C1 C2];


space=linspace(low_point,up_point,tot_point);
      for k=1:tot_point
          G1=C*((1j*space(k)*E-A)\B);
          G2=Cr*((1j*space(k)*Er-Ar)\Br)+Dar;
          G3=Cr_fl*((1j*space(k)*Er_fl-Ar_fl)\Br_fl)+Dar;

          tf(k)=max(svds(G1));
          tf_rom(k)=max(svds(G2));
          tf_rom_fl(k) = max(svds(G3));
          abs_err(k)=max(svds(G1-G2));
          abs_err_fl(k)=max(svds(G1-G3));
          
      end
          rel_err=abs_err./tf;
          rel_err_fl=abs_err_fl./tf;
%  Transfer Function Plot for full,frequency limited and frequency unrestricted reduced model    
        figure(1);
        semilogy(space,tf,'k')
        hold on;
        semilogy(space,tf_rom,'b--')
        hold on;
        semilogy(space,tf_rom_fl,'r-.')
        xlabel('\omega')
        ylabel('\sigma_{max}(G(j\omega))')
        title('Transfer function')
        legend('Full Model','Frequency Unrestricted', 'Frequency Restricted')
        hold off
%  Absolute Error
        figure(2);
       semilogy(space,abs_err,'b--')
        hold on;
        semilogy(space, abs_err_fl, 'r-.')
         xlabel('\omega')
         ylabel('\sigma_{max}(G(j\omega)-G_r(j\omega))')
         title('absolute error')
         legend('Frequency unrestricted', 'Frequency Restricted')
         hold off
%   relative error;
       figure(3);
       semilogy(space, rel_err,'b--')
       hold on;
       semilogy(space,rel_err_fl,'r-.')
       xlabel('\omega')
       ylabel('\sigma_{max}(G(j\omega)-G_r(j\omega))/ \sigma_{max}(G(j\omega))')
       title('Relative error')
       legend('Frequency Unrestricted', 'Frequency Restricted')
 hold off
end