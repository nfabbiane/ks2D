clc, clear all, close all
addpath('../')

%% Parameters

% Kuramoto
omegarmax = 2.5e-3; % maximum growth rate
alphamax = 1.5e-1;  % alpha of maximum growth rate
betamax = 1.5e-1;   % maximum unstable beta

V  = 1.0;           % phase speed (= group speed, KS is not dispersive)

P = 2 * alphamax^2;
R = P^2/(4*omegarmax);
S  = omegarmax*R/betamax^4;

% space discretisation
NX =  72;           % number of modes in x
NZ =  12;           % number of modes in z
LX = 500;           % domain length (x)
LZ = 200;           % domain width (z)
Lf = 150;           % fringe length (x)




%% Initialization
[A,xx,zz,~,alpha,beta] = ks_init(P,R,S,V,LX,LZ,Lf,NX,NZ);




%% Test

% generate field
[q,f,v] = ks_init_input([0 0],[10 10],xx,zz);

% generate output matrix
[C,Cfou,Cphy] = ks_init_output([0 0],[10 10],xx,zz,LX,LZ);

% compute output
dx = xx(2,1) - xx(1,1);
dz = zz(1,2) - zz(1,1);

y    = C*q
yfou = sum(sum(Cfou .* f)) * (LX*LZ)
yphy = sum(sum(Cphy .* v)) * (dx*dz)