clear ; 
clc; 
close all
%% This script test frequency-limited 
%% H2 optimal model reduction generalized system
%%
%load('/home/monir/scratch/software/data/1st_order/dae/index2/NSE_unstable/NSE_unstable_1087.mat');
load('/home/monir/scratch/software/data/1st_order/dae/index2/stokes/stokes.mat')
I=speye(size(A1,1));
PI=I-A2*((A2'*(E1\A2))\(A2'/E1));
%A=PI*A1*PI';E=PI*E1*PI';B=PI*B1;C=C1*PI';
%[U,S,V]=svd(full(PI));
%Ql=U(:,1:825)/diag(sqrt(diag(S(1:825,1:825))));
%Qr=diag(sqrt(diag(S(1:825,1:825))))\V(:,1:825)';
%At=Qr*A1*Qr';Et=Qr*E1*Qr';Bt=Qr*B1;Ct=C1*Qr';
At=PI*A1;Et=E1;Bt=PI*B1;Ct=C1;
%load('/home/monir/scratch/software/data/1st_order/generalized/tchain.mat')%%data_gen.mat'
%%
% load 1storder_butterflygyro.mat
%E=eye(size(A,1));
r = 50; 
maxiter = 50; 
tol = 1e-10;
%% IRKA for initial ROM
%--------------------------------
%[Ar, Br, Cr] = IRKA(Et,At,Bt,Ct,r,maxiter,tol);
[Ar,Br,Cr] = IRKA_index2_nonsymmetric(E1,A1,A2,B1,C1,r,maxiter,tol);
%% TSIA for frequency unrestricted
%----------------------------------------
%[Ar,Br,Cr] = gen_TSIA(Et,At,Bt,Ct,Ar,Br,Cr,maxiter,tol);
[Ar1,Br1,Cr1] = index2_TSIA(E1,A1,A2,B1,C1,Ar,Br,Cr,maxiter,tol);
%% Frequency range [w1,w2]
%-------------------------------
---
w1=2; 
w2=4;
[Ar2, Br2, Cr2] = gen_Fl_TSIA(E1,A1,A2,B1,C1,Ar,Br,Cr,maxiter,tol,w1,w2);
%[Ar2,Br2,Cr2] = index2_Fl_TSIA(E1,A1,A2,B1,C1,Ar,Br,Cr,maxiter,tol,w1,w2);
%%
low_point =1;      %Lower point of the plotting domain
up_point = 6;       % Higher point of the plotting domain
tot_point = 500;   
Er= speye(size(Ar,1));
Er_fl= Er;
%% Ploat transfer function
%--------------------------------
[space,tf] = tf_plot(At,Bt,Ct,Et,Ar1,Br1,Cr1,Er,Ar2,Br2,Cr2,Er,low_point,up_point,tot_point);
%%
%[space,tf] =sigma_plot(A,B,C,E,E,A,B,C,low_point,up_point,tot_point);