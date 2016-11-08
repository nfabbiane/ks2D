clc, clear all, close all
addpath('../lib','../')

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
[A,xx,zz] = ks_init(P,R,S,V,LX,LZ,Lf,NX,NZ);




%% Test

% generate field
[q,f,v] = ks_init_input([0 0],[20 20],xx,zz);

% use q2v and v2q
v1 = q2v(q,NX,NZ);
v2 = q2v(v2q(v),NX,NZ);

% plot result
figure(1001); clf       
subplot(2,1,1); surf(xx,zz,v1-v,'EdgeColor','none');
                colorbar('EO'); colormap(redblue)
                cax = caxis; caxis([-1 1]*max(abs(cax)));
                axis image; view(2);
                xlabel('x'), ylabel('z');
                title('v - q2v(q)')
subplot(2,1,2); surf(xx,zz,v2-v,'EdgeColor','none');
                colorbar('EO'); colormap(redblue)
                cax = caxis; caxis([-1 1]*max(abs(cax)));
                axis image; view(2);
                xlabel('x'), ylabel('z');
                title('v - q2v(v2q(v))')