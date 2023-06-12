clear ; 
clc; 
close all
%% This script test frequency-limited 
%% H2 optimal model reduction generalized system
%%
load('D:/data/1st_order/generalized/tchain.mat')%%data_gen.mat'
%%
%A=E\A;
%B=E\B;
%E=eye(size(A,1));
%load('/home/monir/scratch/software/data/1st_order/dae/index2/NSE_unstable/NSE_unstable_1087.mat');
%I=speye(size(A1,1));
%PI=I-A2*((A2'*(E1\A2))\(A2'/E1));
%A=PI*A1*PI';E=PI*E1*PI';B=PI*B1;C=C1*PI';
%[U,S,V]=svd(full(PI));
%Ql=U(:,1:825)/diag(sqrt(diag(S(1:825,1:825))));
%Qr=diag(sqrt(diag(S(1:825,1:825))))\V(:,1:825)';
%At=Qr*A1*Qr';Et=Qr*E1*Qr';Bt=Qr*B1;Ct=C1*Qr';
%At=PI*A1;Et=E1;Bt=PI*B1;Ct=C1;
%load('/home/monir/scratch/software/data/1st_order/generalized/tchain.mat')%%data_gen.mat'
%%
% load 1storder_butterflygyro.mat
E=speye(size(A,1));
r = 50; 
maxiter = 50; 
tol = 1e-10;
%% IRKA for initial ROM
%--------------------------------
[Ar,Br,Cr,So,bo,co,V,W,err_hist] = IRKA(E,A,B,C,r,maxiter,tol);
%% TSIA for frequency unrestricted
%----------------------------------------
% [Ar,Br,Cr,V,W] = gen_TSIA(E,A,B,C,Ar,Br,Cr,maxiter,tol);
%% Frequency range [w1,w2]
%----------------------------------
w1=1; 
w2=2;
%tic;[Ar_fl, Br_fl, Cr_fl] = gen_Fl_TSIA(E,A,B,C,Ar,Br,Cr,maxiter,tol,w1,w2);toc
tic;[Ar_fl, Br_fl, Cr_fl] = gen_Fl_TSIA_ex(E,A,B,C,Ar,Br,Cr,V,W,maxiter,tol,w1,w2);toc
%%
low_point =0;      %Lower point of the plotting domain
up_point = 4;       % Higher point of the plotting domain
tot_point = 500;   
Er= speye(size(Ar,1));
Er_fl= Er;
%% Ploat transfer function
%--------------------------------
[space,tf] = tf_plot(A,B,C,E,Ar,Br,Cr,Er,Ar_fl,Br_fl,Cr_fl,Er_fl,low_point,up_point,tot_point);
%%
%[space,tf] =sigma_plot(A,B,C,E,E,A,B,C,low_point,up_point,tot_point);