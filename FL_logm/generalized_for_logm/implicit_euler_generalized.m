function [t,yout,youtr,abserr,relerr]=implicit_euler_generalized(E,A,B,C,Er,Ar,Br,Cr,Ert,Art,Brt,Crt,t0,tf,h)
n = abs((tf-t0)/h); % Calculate steps
 t = [t0 zeros(1,n-1)];
 yout = zeros(size(C,1),n);
 youtr = zeros(size(Cr,1),n);
 [L,U,P,Q] = lu(E-h*A);
% [Lr,Ur,Pr,Qr] = lu(Er-h*Ar);
 x=zeros(size(A,1),1);
 xr=zeros(size(Ar,1),1);
  xrt=zeros(size(Ar,1),1);
 for i=1:n-1
      t(i+1) = t(i)+h;
%       if i<=20
%             x=Q*(U\(L\(P*(E*x))));
%             xr=(Er-h*Ar)\(Er*xr);
%              xrt=(Ert-h*Art)\(Ert*xr);
%       else
    % x=Q*(U\(L\(P*(E*x+h*sum(B,2)))));
    % xr=(Er-h*Ar)\(Er*xr+h*sum(Br,2));
    x=Q*(U\(L\(P*(E*x+h*B))));
    xr=(Er-h*Ar)\(Er*xr+h*Br);
    xrt=(Ert-h*Art)\(Ert*xrt+h*Brt);
    % end
     yout(:,i+1)=C*x;
     youtr(:,i+1)=Cr*xr;
     youtrt(:,i+1)=Crt*xrt;
 end
 abserr=abs(yout-youtr);
  abserrt=abs(yout-youtrt);
 relerr=abserr./abs(yout);
  relerrt=abserrt./abs(yout);
     figure(200)
     plot(t,yout)
    % legend('original model')
     hold on 
     plot(t,youtr,'r-.')
     plot(t,youtrt,'m-.')
     legend('original','stnd-irna','tl-irka')
    title('TF')
     hold off
    figure(205)
    semilogy(t,abserr)
    hold on
    semilogy(t,abserrt,'r-.')
    legend('stnd-irna','tl-irka')
    title('absolute error')
    hold off
    figure(206)
    semilogy(t,relerr)
    hold on
    semilogy(t,relerrt,'r-.')
    title('relative error')
    legend('stnd-irna','tl-irka')
    hold off
end