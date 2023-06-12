clear;
clc; 
close all
%%
load('/home/monir/scratch/software/data/1st_order/dae/index1/bps_1142.mat')
%%
Dar = -C2*(J4\B2);
%%
E1 = speye(size(J1,1));
%% 
r = 40; 
maxiter = 110; 
tol = 1e-8;
%%
tic;[Ar,Br,Cr] = IRKA_index1(E1,J1,J2,J3,J4,B1,B2,C1,C2,r,maxiter,tol);toc
%%
%[Ar,Br,Cr] = index1_TSIA(E1,J1,J2,J3,J4,B1,B2,C1,C2,Ar,Br,Cr,maxiter,tol);
%%
w1=4; w2=6;
%% FL-TSIA
%[Ar_fl,Br_fl,Cr_fl,S] = IRKA_fl_index1(E1,J1,J2,J3,J4,B1,B2,C1,C2,Ar, Br, Cr,r,maxiter,tol,w1,w2);
[Ar_fl,Br_fl,Cr_fl] = index1_Fl_TSIA(E1,J1,J2,J3,J4,B1,B2,C1,C2,Ar,Br,Cr,maxiter,tol,w1,w2);
%%
low =2;      % Lower point of the plotting domain
up = 10;     % Higher point of the plotting domain
points = 200;
%%
Er=speye(size(Ar,1));
Er_fl= Er;
[space,tf] = tf_plot(E1,J1,J2,J3,J4,B1,B2,C1,C2,Er,Ar,Br,Cr,Er_fl,Ar_fl,Br_fl,Cr_fl,Dar,low,up,points);
