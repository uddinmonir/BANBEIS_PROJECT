function [space,tf] =tf_plot(A,B,C,E,Ar,Br,Cr,Er,Ar_fl,Br_fl,Cr_fl,Er_fl,low_point,up_point,tot_point)
     space=linspace(low_point,up_point,tot_point);
     
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
        figure(10);
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
        figure(20);
       semilogy(space,abs_err,'b--')
        hold on;
        semilogy(space, abs_err_fl, 'r-.')
         xlabel('\omega')
         ylabel('\sigma_{max}(G(j\omega)-G_r(j\omega))')
         title('absolute error')
         legend('Frequency unrestricted', 'Frequency Restricted')
         hold off
%   relative error;
       figure(30);
       semilogy(space, rel_err,'b--')
       hold on;
       semilogy(space,rel_err_fl,'r-.')
       xlabel('\omega')
       ylabel('\sigma_{max}(G(j\omega)-G_r(j\omega))/ \sigma_{max}(G(j\omega))')
       title('Relative error')
       legend('Frequency Unrestricted', 'Frequency Restricted')
 hold off
end