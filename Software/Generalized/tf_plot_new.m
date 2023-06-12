function [space,tf] =tf_plot_new(A,B,C,E,Ar,Br,Cr,Er,low_point,up_point,tot_point)
     space=linspace(low_point,up_point,tot_point);
     
      for k=1:tot_point
          G1=C*((1j*space(k)*E-A)\B);
          G2=Cr*((1j*space(k)*Er-Ar)\Br);

          tf(k)=max(svds(G1));
          tf_rom(k)=max(svds(G2));
          abs_err(k)=max(svds(G1-G2));
          
      end
          rel_err=abs_err./tf;
%  Transfer Function Plot for full,frequency limited and frequency unrestricted reduced model    
        figure(1);
        semilogy(space,tf,'r')
        hold on;
        semilogy(space,tf_rom,'b--')
        xlabel('\omega')
        ylabel('\sigma_{max}(G(j\omega))')
        title('Transfer function')
        legend('Full Model','Frequency Unrestricted', 'Frequency Restricted')
        hold off
%  Absolute Error
        figure(2);
       semilogy(space,abs_err,'b--')
         xlabel('\omega')
         ylabel('\sigma_{max}(G(j\omega)-G_r(j\omega))')
         title('absolute error')
         legend('Frequency unrestricted', 'Frequency Restricted')
         hold off
%   relative error;
       figure(3);
       semilogy(space, rel_err,'b--')
       xlabel('\omega')
       ylabel('\sigma_{max}(G(j\omega)-G_r(j\omega))/ \sigma_{max}(G(j\omega))')
       title('Relative error')
       legend('Frequency Unrestricted', 'Frequency Restricted')
 hold off
end