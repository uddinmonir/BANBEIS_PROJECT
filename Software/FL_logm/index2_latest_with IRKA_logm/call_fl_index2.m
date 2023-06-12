clear all;
clc;
clf;
close all
%% This script test frequency-limited
%% H2 optimal model reduction generalized system
%%
%load('NSE1.mat')
Dim = 3;
data_NSE_n;

I=speye(size(A1,1));
PI=I-A2*((A2'*(E1\A2))\(A2'/E1));
At=PI*A1;Et=E1;Bt=PI*B1;Ct=C1;
r = 120;
maxiter1 = 50;
maxiter2 = 50;
tol = 1e-10;
%% IRKA for initial ROM
%--------------------------------
%[Ar, Br, Cr] = IRKA(Et,At,Bt,Ct,r,maxiter,tol);
tic;
[Ar1,Br1,Cr1,V1] = IRKA_index2_nonsymmetric(E1,A1,A2,B1,C1,r,maxiter1,tol);
%% TSIA for frequency unrestricted
%----------------------------------------
%[Ar,Br,Cr] = gen_TSIA(Et,At,Bt,Ct,Ar,Br,Cr,maxiter,tol);
%[Ar2,Br2,Cr2,V,W] = index2_TSIA(E1,A1,A2,B1,C1,Ar1,Br1,Cr1,maxiter1,tol);
%% Frequency range [w1,w2]
%----------------------------------
r1=500;
w1=4;
w2=5;
%[Ar2, Br2, Cr2] = index2_Fl_TSIA_FB(K1,L1,E1,A1,A2,B1,C1,Ar1,Br1,Cr1,r1,maxiter,tol,w1,w2);
[Ar2, Br2, Cr2] = index2_Fl_TSIA(E1,A1,A2,B1,C1,Ar1,Br1,Cr1,r1,maxiter2,tol,w1,w2);
%[Ar2,Br2,Cr2] = index2_Fl_irka(E1,A1,A2,B1,C1,Ar1,Br1,Cr1,r1,maxiter,tol,w1,w2);
%%
toc;
low_point=1;      %Lower point of the plotting domain
up_point=6;       % Higher point of the plotting domain
tot_point=100;
Er1=speye(size(Ar1,1));
Er2=Er1;
%% Ploat transfer function
tic;[space,tf] =tf_plotex(E1,A1,A2,B1,C1,Er1,Ar1,Br1,Cr1,Er2,Ar2,Br2,Cr2,low_point,up_point,tot_point);toc
%--------------------------------
