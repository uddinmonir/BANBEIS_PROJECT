function [Vn]=rksm(E,A,B,s,m,tol,tolY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reformulation
%A=J1-J2*(J4\J3);
[n,~]=size(A);
%E=speye(n); %b=B1-J2*(J4\B2);%B=full(b);
p=size(B,2);
V=(A-s*E)\B; %V=V(1:n,:);  
[V,~]=qr(V,0);
Br=V'*B;
I=speye(p);
O=0*I;  %Vn=zeros(n,p*(m+2));%Vn(1:n,1:p)=V;
Vn=V;
H=zeros(p*(m+2),p*(m+1));
resnorm=[];
%%
nrma=norm(A,'fro');
[~,w]=qr(B,0);
ww=inv(w);
nrmb=norm(inv(ww),'fro')^2; 
%%
Ar=V'*A*V; %V'AVr  %Br=(V'*B)*(B'*V);
i=1;
j=0;
Bn=sprand(n,100,.8);
[Q,~]=qr(full(Bn),0);
sp=adaptive_ADIshift(E,A,Q,5);  %sp=eigs(Ar,4); % compute l ADI shift parameters
l=length(sp);
%%%%%%%%%%%%%%%%%%%%%%%%%
while i < m
  %i=i+1;
     if j<l
          j=j+1;
       else              % k=(i*p-l*p)+1;         %Vnn=Vn(1:n,k:i*p);
         m2=size(Vn,2);
         Vcn=Vn(:,(m2-l*p)+1:end);    %size(Vcn,2)
         [Q,~]=qr(full(Vcn),0);          % sp=eigs(Ar,16);%%% next shifts computation
         sp = adaptive_ADIshift(E,A,Q,4);
         l=length(sp);
         j=1;
     end
   s=sp(j); 
   if real(s)<0,
       s=-real(s)+sqrt(-1)*imag(s);
   end
 %%   
    Atil=(A-s*E);      %Vtil=[V;sparse(size(J4,1),size(V,2))];
    wrk = Atil\V;     % wrk= wrk(1:n,:);  
%% Gram-Schmidt step
 jms=(i-1)*p+1;  j1s=(i+1)*p;   js=i*p;  js1=js+1;
  for it=1:2,
      for k=1:i
        k1=(k-1)*p+1;
        k2=k*p; 
        gamma=Vn(1:n,k1:k2)'*wrk;
        H(k1:k2,jms:js)=H(k1:k2,jms:js)+gamma;
        wrk=wrk-Vn(:,k1:k2)*gamma;
      end
  end
    [V,H(js1:j1s,jms:js)]=qr(wrk,0);
   %% Solve Lyapunov equation  % b=speye(js,p)*B*speye(js,p)';
  Y=lyap(Ar,(Br*Br'));  
  nrmx=norm(Y,'fro');
   newAv=A*V;
   g = Vn(1:n,1:js)'*newAv;
   u1=newAv-Vn(1:n,1:js)*g;
   d=-Vn(1:n,1:js)*(Y*(H(1:i*p,1:i*p)'\[sparse(p*(i-1),p);I])*H(js1:j1s,jms:js)');
   U=[-V*s(end), d u1 ];
   rr=qr(full(U),0); 
   rr=triu(rr(1:size(rr,2),:));
   k=i;
   normres(k)=norm(rr*sparse([O I O; I O I; O I O ])*rr','fro')/(nrmb+nrma*nrmx);
   resnorm=[resnorm,normres];
    fprintf(1,'step: %4d norm. resid.: %d\n',k,normres(k))
   if(normres(k)<tol)   
     break
   end 
     g1=g; 
     g2=V'*A*Vn(1:n,1:js);
     g3=V'*A*V;
     Ar=[Ar g1; g2, g3];
     Bn=V'*B;
     Br=[Br;Bn]; % size(Ar)%  size(Br) % Vn(1:n,js1:j1s)=V; 
    Vn=[Vn,V];   % size(Vn,2)
   i=i+1;
 end;
%% Reduce rank of solution, if needed
[uY,sY]=eig(Y); 
[sY,id]=sort(diag(sY));
sY=flipud(sY);
uY=uY(:,id(end:-1:1));
is=sum(abs(sY)>tolY);
Y0 = uY(:,1:is)*diag(sqrt(sY(1:is))); 
Z = Vn(:,1:size(Y0,1))*Y0;
Final_Z=size(Z)
RKStotal_time=toc;
fprintf('Space dim %d  Solution rank %d  time %d \n',(i+1)*p,is,RKStotal_time);
return