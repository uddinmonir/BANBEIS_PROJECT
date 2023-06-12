clear ;
clc;
clf;
close all;
%% This script test frequency-limited
%% H2 optimal model reduction generalized system
%%
%load('/home/monir/scratch/software/data/1st_order/generalized/iss.mat')%%data_gen.mat'
%%
%load 1storder_butterflygyro.mat
%load ('D:/data/1st_order/generalized/tchain.mat')
load('D:/data/1st_order/dae/index1/NSE_unstable_1087.mat');
PI = speye(size(A1,1)) - A2*((A2'*(M\A2))\(A2'/M));
%E=eye(size(A,1));

EP = PI*M*PI';
AP = PI*A1*PI';
BP = PI*B1;
CP = C1*PI';
r  = 50;
maxiter = 50;
tol = 1e-8;

low_point =0;      %Lower point of the plotting domain
up_point = 6;      % Higher point of the plotting domain
tot_point = 1000;
% Er= speye(size(Ar,1));
% Er_fl= Er;
[space,tf] = tf_plot_new(A,B,C,E,AP,BP,CP,EP,low_point,up_point,tot_point);
%% IRKA for initial ROM
%--------------------------------
% [Ar, Br, Cr] = IRKA(E,A,B,C,r,maxiter,tol);
% %% TSIA for frequency unrestricted
% %----------------------------------------
% [Ar,Br,Cr] = gen_TSIA(E,A,B,C,Ar,Br,Cr,maxiter,tol);
% %% Frequency range [w1,w2]
% %----------------------------------
% w1=2;
% w2=4;
% tic;[Ar_fl, Br_fl, Cr_fl] = gen_Fl_TSIA(E,A,B,C,Ar,Br,Cr,maxiter,tol,w1,w2);toc
% %%
% low_point =0;      %Lower point of the plotting domain
% up_point = 6;       % Higher point of the plotting domain
% tot_point = 1000;
% Er= speye(size(Ar,1));
% Er_fl= Er;
%% Ploat transfer function
%--------------------------------
%[space,tf] = tf_plot(A,B,C,E,Ar,Br,Cr,Er,Ar_fl,Br_fl,Cr_fl,Er_fl,low_point,up_point,tot_point);
%%
%[space,tf] =sigma_plot(A,B,C,E,Ar,Br,Cr,Er,Ar_fl,Br_fl,Cr_fl,Er_fl,low_point,up_point,tot_point);
% E_e = blkdiag(E,Er);
% E_fl = blkdiag(E,Er_fl);
% %
% A_e = blkdiag(A,Ar);
% A_fl = blkdiag(A,Ar_fl);
% %
% B_e = [B;Br];  B_fl = [B;Br_fl];
% C_e = [C -Cr]; C_fl = [C -Cr_fl];
% %
% S_omega= real((1i/pi)*logm(full((A+1i*w1*E)\(A+1i*w2*E))));
% S_omega_e= real((1i/pi)*logm(full((A_e+1i*w1*E_e)\(A_e+1i*w2*E_e))));
% S_omega_fl= real((1i/pi)*logm(full((A_fl+1i*w1*E_fl)\(A_fl+1i*w2*E_fl))));
% %S_omega = real((1i/pi)*logm(full(S_omega)));
% P = S_omega*(B*B')+(B*B')*S_omega';
% P_e = S_omega_e*(B_e*B_e')+(B_e*B_e')*S_omega_e';
% P_fl = S_omega_fl*(B_fl*B_fl')+(B_fl*B_fl')*S_omega_fl';
% Q = S_omega'*(C'*C)+(C'*C)*S_omega;
% Q_e = S_omega_e'*(C_e'*C_e)+(C_e'*C_e)*S_omega_e;
% Q_fl = S_omega_fl'*(C_fl'*C_fl)+(C_fl'*C_fl)*S_omega_fl;
% %
% FullH2 = sqrt(trace(C*lyap(A,P,[],E)*C'));
% FullH2_IRKA = sqrt(trace(C_e*lyap(A_e,P_e,[],E_e )*C_e'));
% FullH2_IRKA_fl = sqrt(trace(C_fl*lyap(A_fl,P_fl,[],E_fl)*C_fl'));
% %
% ErrH2IRKA = (FullH2_IRKA)/FullH2
% ErrH2IRKA_fl = (FullH2_IRKA_fl)/FullH2




