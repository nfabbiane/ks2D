clc, close all, clear all

%% Parameters

% Kuramoto
alphamax = 1.68e-1;
omegarmax = 2.67e-3;
betamax = 2.15e-1;

P = 2 * alphamax^2;
R = P^2/(4*omegarmax);
S  = omegarmax*R/betamax^4;
V  = 0.4;

% space discretisation
nqx =   96;
nqz =   16;
Lx  = 1000;
Lz  =  150;
Lf  =  200;


%% Grid
alpha0 = 2*pi / Lx;
ialpha = -floor(nqx/2):ceil(nqx/2-1)';
alpha  = alpha0 * ialpha;

beta0 = 2*pi / Lz;
ibeta = -floor(nqz/2):ceil(nqz/2-1)';
beta  = beta0 * ibeta;

[beta,alpha] = meshgrid(beta,alpha);


%% Dispersion relation
omegar = (P*(alpha.^2) - (alpha.^4 + S*beta.^4))/R;
omegai = V*alpha.*(1 + beta.^2/(8*P));


%% Plot
ax = [[-1 1]*1.5*betamax [0 1]*2*alphamax];
figure(1); clf;

subplot(2,1,1); contourf(beta(1,:),alpha(:,1),omegar,(-1:.01:1)*20*omegarmax,'LineColor','none'); hold on
                contour(beta(1,:),alpha(:,1),omegar,[0 0],'-w','Linewidth',2); hold off
                axis image; axis(ax); colorbar('EO');
                cax = caxis; caxis([-1 1]*max(cax));
                xlabel('\beta'); ylabel('\alpha'); title('\omega_r')

subplot(2,1,2); contourf(beta(1,:),alpha(:,1),omegai,(0:.1:1)*V*alphamax,'LineColor','none'); 
                axis image; colorbar('EO');
                cax = caxis; axis(ax); caxis([0 1]*max(abs(cax)));
                xlabel('\beta'); ylabel('\alpha'); title('\omega_r')