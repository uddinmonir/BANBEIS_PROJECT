function p = adaptive_ADIshift(E,A,V,l)
% l:= desirded number of  ADI-shifts
% function Shift_Parameter_new compute desired number of 
% ADI shift parameters for the system 
% 
%       Ex'(t) = Ax(t)+Bu(t) (Index-2 system)
%
%  where A=[ A11 A12 ]  E = [ M  0 ]  B= [ B1 ]
%          [ A21^T 0 ],     [ 0  0 ],    [ 0  ]
% 
%     MOHAMMAD MONIR UDDIN (May 2013)
%  n1=size(A,1);
%  B=sprand(n1,100,.8);
%  [Q,~]=svd(full(B),0);
%   s=diag(S);
%  tol = max(size(A)) * eps(max(s));
%  r = sum(s > tol);
%  Q=Uf(:,1:r)*S(1:r,1:r);
   En=(V'*E*V);
   An=(V'*A*V);
   rw=eig(full(An));
    rw0 = rw;
    rw = [];
   for j = 1:length(rw0)
       if real(rw0(j))<0
         rw=[rw;rw0(j)];
       else 
         rw=[rw;-real(rw0(j))+sqrt(-1)*imag(rw0(j))];
       end
   end
   p=lp_mnmx(rw,l); % solve the minimax problem
   p=-real(p)+sqrt(-1)*imag(p);
    %p=rw; % solve the minimax problem
    % p=-real(p)+sqrt(-1)*imag(p);
  % p=rw(1:l,:);
end