function [t,yout,youtr,abserr,relerr]=implicit_euler_indexI(E,A,B,C,D,ntEr,ntAr,ntBr,ntCr,Ert,Art,Brt,Crt,t0,tf,h)
n = abs((tf-t0)/h); % Calculate steps
 t = [t0 zeros(1,n-1)];
 yout = zeros(size(C,1),n);
 youtr = zeros(size(ntCr,1),n);
 youtrt = zeros(size(Crt,1),n);
 [L,U,P,Q] = lu(E-h*A);
 x=zeros(size(A,1),1);
 xr=zeros(size(ntAr,1),1);
  xrt=zeros(size(Art,1),1);
 for i=1:n-1
      t(i+1) = t(i)+h;
    x=Q*(U\(L\(P*(E*x+h*B))));
    xr=(ntEr-h*ntAr)\(ntEr*xr+h*ntBr);
    xrt=(Ert-h*Art)\(Ert*xrt+h*Brt);
     yout(:,i+1)=C*x;
     youtr(:,i+1)=ntCr*xr;
     youtrt(:,i+1)=Crt*xrt;
 end
 abserr=abs(yout-youtr);
  abserrt=abs(yout-youtrt);
 relerr=abserr./abs(yout);
  relerrt=abserrt./abs(yout);
     figure(1)
     plot(t,yout)
     hold on 
     plot(t,youtr,'r--')
     plot(t,youtrt,'m-.')
     legend('original','Time Unlimited','Time Limited')
    title('TF')
     hold off
    figure(2)
    semilogy(t,abserr,'r--')
    hold on
    semilogy(t,abserrt,'m-.')
    legend('Time Unlimited','Time Limited')
    title('absolute error')
    hold off
    figure(3)
    semilogy(t,relerr,'r--')
    hold on
    semilogy(t,relerrt,'m-.')
    title('relative error')
    legend('Time Unlimited','Time Limited')
    hold off
end