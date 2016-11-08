clc, clear all, %close all
addpath('./lib')

%% Parameters

% Kuramoto
omegarmax = 2.5e-3; % maximum growth rate
alphamax = 1.5e-1;  % alpha of maximum growth rate
betamax = 1.5e-1;   % maximum unstable beta

V  = 1.0;           % phase speed (= group speed, KS is not dispersive)

P = 2 * alphamax^2;
R = P^2/(4*omegarmax);
S = omegarmax*R/betamax^4;

% space discretisation
NX =  72;           % number of modes in x
NZ =  12;           % number of modes in z
LX = 500;           % domain length (x)
LZ = 200;           % domain width (z)
Lf = 150;           % fringe length (x)




%% Initialization
[~,~,~,omega,alpha,beta] = ks_init(P,R,S,V,LX,LZ,Lf,NX,NZ);

% - reorder for better contours
alpha = fftshift(alpha);
beta  = fftshift(beta );
omega = fftshift(omega);




%% Plot
lvlr = (-1:.1:1) * omegarmax;
lvli = (-1:.1:1) * V * 2*alphamax;

figure(1); clf;

subplot(2,1,1); contourf(alpha(:,1),beta(1,:),real(omega).',lvlr,'LineColor','none');
                axis image; colorbar('EO'); colormap(gca,redblue);
                xlabel('\alpha'); ylabel('\beta'); title('\omega_r')

subplot(2,1,2); contourf(alpha(:,1),beta(1,:),imag(omega).',lvli,'LineColor','none'); 
                axis image; colorbar('EO');
                xlabel('\alpha'); ylabel('\beta'); title('\omega_i')