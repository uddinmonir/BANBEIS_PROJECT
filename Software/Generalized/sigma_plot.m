function [space,tf] =sigma_plot(A,B,C,E,Ar,Br,Cr,Er,Ar_fl,Br_fl,Cr_fl,Er_fl,low_point,up_point,tot_point)
     space=logspace(low_point,up_point,tot_point);
     
      for k=1:tot_point
          G1=C*((1j*space(k)*E-A)\B);
          G2=Cr*((1j*space(k)*Er-Ar)\Br);
          G3=Cr_fl*((1j*space(k)*Er_fl-Ar_fl)\Br_fl);

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
        loglog(space,tf,'k')
        hold on;
        loglog(space,tf_rom,'b--')
        hold on;
        loglog(space,tf_rom_fl,'r-.')
        xlabel('\omega')
        ylabel('\sigma_{max}(G(j\omega))')
        title('Transfer function')
        legend('Full Model','Frequency Unrestricted', 'Frequency Restricted')
        hold off
%  Absolute Error
        figure(2);
       loglog(space,abs_err,'b--')
        hold on;
        loglog(space, abs_err_fl, 'r-.')
         xlabel('\omega')
         ylabel('\sigma_{max}(G(j\omega)-G_r(j\omega))')
         title('absolute error')
         legend('Frequency unrestricted', 'Frequency Restricted')
         hold off
%   relative error;
       figure(3);
       loglog(space, rel_err,'b--')
       hold on;
       loglog(space,rel_err_fl,'r-.')
       xlabel('\omega')
       ylabel('\sigma_{max}(G(j\omega)-G_r(j\omega))/ \sigma_{max}(G(j\omega))')
       title('Relative error')
       legend('Frequency Unrestricted', 'Frequency Restricted')
 hold off
end