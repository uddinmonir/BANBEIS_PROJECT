% % -----------------------------------------------------
% % stable data::::of Navier-Stokes Equation::::::
% %-----------------------------------------------
% %load nse_cell1.mat
% load re500_all.mat %( old data)
 load mat_nse_re_500_new.mat % (new data)
% %  Ans: (write mat.mat_v to find follows)
% %  E: {1x6 cell}
% %          M: {1x6 cell}
% %      fullA: {1x6 cell}
% %          A: {1x6 cell}
% %          S: {1x6 cell}
% %          K: {1x6 cell}
% %          R: {1x6 cell}
% %          G: {1x6 cell}
% %         Re: 500
% %          B: {[3452x2 double] [8726x2 double]  [20512x2 double]  
%               [45718x2 double] [99652x2 double]  [211452x2 double]}
% %          C: {[7x3452 double]  [11x8726 double]  [17x20512 double]  
%                [23x45718 double]  [37x99652 double]  [47x211452 double]}
% % %     Feed_0: {[3452x2 double]  [8726x2 double]  [20512x2 double]}
 % Dim=5;
 A1=mat.mat_v.A{Dim};
 A2=mat.mat_v.G{Dim};
 [n1,n2]=size(A2);
 A3=A2';
 A4=sparse(n2,n2);
 E1=mat.mat_v.M{Dim};
 B1=mat.mat_v.B{Dim};
 C1=mat.mat_v.C{Dim};
 K1=real(mat.mat_v.Feed_0{Dim});
 B2=sparse(n2,size(B1,2));
 C2=sparse(size(C1,1),n2);
 K2=B2;
%  A_til=[J1_til A2;A2' A4]-[B1;B2]*[k;B2]';
%  Atil=A1-B1*K1';
%  A1=A1-.5*M; %shift eigenvalues from imaginary to left
   E=[E1 sparse(n1,n2);sparse(n2,n1) A4];
   A=[A1 A2;A3 A4];
% Atil=[A1-B1*K1' A2;A2' A4];
% Ast=[A1 A2;A2' A4];
  B=[B1;B2];
  C=[C1 C2];
  K=[K1; K2];
% Da=spalloc(size(C2,1),size(B1,2),0);
% I=speye(size(A1,1));
% PI=I-A2*((A2'*(M\A2))\(A2'*(M)^(-1)));
% A1=J1-J2*((J3*(E1\J2))\(J3*(E1\J1)));