clear ; clc; close all

%% Define random matrice

E=E1;
GJ1 = E\J1;
GJ2 = E\J2;
GB1 = E\B1;
JJ43 = J4\J3; 
JB42 = J4\B2;
A = GJ1 - (GJ2*JJ43);
B = GB1 - (GJ2*JB42);
C = C1 - (C2*JJ43);
D = - (C2*JB42);
n = size(A,1);
m = size(B,2);
p = size(C,1);
r = 50; 
Imaxiter = 50; 
Tmaxiter = 50;
TLmaxiter = 50;
tol = 1e-8;
Time1 = 0;
Time2 = 4;
Ctol = 1e-16;    % Convergency limit 
nZtol=1e-20;  % Non zero tolerance value
is = 1; % Input Signal
os = 1; % Output Signal
[Ar, Br, Cr] = IndexIIRKA(GJ1,GJ2,JJ43,GB1,JB42,J3,J4,B2,C1,C2,r,Imaxiter,tol,n,m,p); 

[ntAr,ntBr,ntCr] = indexI_TSIA(A,B,C,Ar,Br,Cr,GJ1,GJ2,JJ43,JB42,GB1,C1,C2,J3,J4,Tmaxiter,tol);

[Art,Brt,Crt] = indexI_tl_TSIA(A,B,C,Ar,Br,Cr,GJ1,GJ2,JJ43,JB42,GB1,C1,C2,J3,J4,TLmaxiter,tol,Time1,Time2,Ctol,nZtol);
t0 = 0;
tf = 6;
h=0.01;
[t,yout,youtr,abserr,relerr]=implicit_euler_indexI(E,A,B(:,is),C(os,:),D,eye(size(ntAr,1)),ntAr,ntBr(:,is),ntCr(os,:),eye(size(Art,1)),Art,Brt(:,is),Crt(os,:),t0,tf,h);
