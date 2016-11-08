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

% time integration
tend = 2500;        % final time
dt   = .25;         % time-step

tscr = tend/20;     % time interval for screen output




%% Initialization
[A,xx,zz,L,alpha,beta,lambda,LAMBDA] = ks_init(P,R,S,V,LX,LZ,Lf,NX,NZ);




%% Inputs matrix B

% disturbance d (Gaussian shape at x_d, z_d with sigma_d variance)
nd = 1; 
posd = zeros(nd,2); posd(:,1) = 0;
                    posd(:,2) = -LZ/2 + LZ/(2*nd):LZ/(nd):LZ/2 - LZ/(2*nd);
sigd = zeros(nd,2); sigd(:,1) = 4;
                    sigd(:,2) = 4;

[~,Bd] = ks_init_input(posd,sigd,xx,zz);

      
% actuator u (Gaussian shape at x_u, z_u with sigma_u variance)
nu = 4; 
posu = zeros(nu,2); posu(:,1) = 200;
                    posu(:,2) = -LZ/2 + LZ/(2*nu):LZ/(nu):LZ/2 - LZ/(2*nu);
sigu = zeros(nu,2); sigu(:,1) = 4;
                    sigu(:,2) = 4;

[~,Bu] = ks_init_input(posu,sigu,xx,zz);




%% Outputs matrix C

% measurement y (Gaussian shape at x_y, z_y with sigma_y variance)
ny = nu; 
posy = zeros(ny,2); posy(:,1) = 100;
                    posy(:,2) = -LZ/2 + LZ/(2*ny):LZ/(ny):LZ/2 - LZ/(2*ny);
sigy = zeros(ny,2); sigy(:,1) = 4;
                    sigy(:,2) = 4;

[~,Cy] = ks_init_output(posy,sigy,xx,zz,LX,LZ);


% output z (Gaussian shape at x_z with sigma_z variance)
nz = nu; 
posz = zeros(nz,2); posz(:,1) = 300;
                    posz(:,2) = -LZ/2 + LZ/(2*nz):LZ/(nz):LZ/2 - LZ/(2*nz);
sigz = zeros(nz,2); sigz(:,1) = 4;
                    sigz(:,2) = 4;

[~,Cz] = ks_init_output(posz,sigz,xx,zz,LX,LZ);




%% Time Integration
t = 0:dt:tend; nt = length(t);

% init state space variables
nq = size(A,1);

q = zeros(NX,NZ,1);
v = zeros(NX,NZ,nt);
f = zeros(NX,NZ,1);

% init signals
d = zeros(nd,nt);
u = zeros(nu,nt);
y = zeros(ny,nt);
z = zeros(nz,nt);

% init disturbaces
% - impulse
d(:,1) = 1/dt;
% - noise
% d(:,:) = randn(nd,nt); % unitary variance Gaussian white noise
% for j = 1:nd
%     d(j,:) = d(j,:) - mean(d(j,:),2);   % enforce zero-mean
%     d(j,:) = d(j,:) / std(d(j,:),[],2); % enforce unitary variance
% end

% time loop
fprintf('\nKS time-integration.\n')

for i = 1:nt-1
    
    % input(s)
    f(:,:) = 0;
    for j = 1:nd % separate for d and u if nd ~= nu
        f = f + Bd(:,:,j) * d(j,i) + Bu(:,:,j) * u(j,i);
    end
    
    
    % KS time-step
    [q(:,:),v(:,:,i+1)] = ks_timestep(q,f,dt,L,lambda);
    
    
    % output(s)
    for j = 1:ny
        y(j,i+1) = real(sum(sum(Cy(:,:,j) .* q(:,:)))) * (LX*LZ);
    end
    
    for j = 1:nz
        z(j,i+1) = real(sum(sum(Cz(:,:,j) .* q(:,:)))) * (LX*LZ);
    end
    
    
    
    
    % simulation status
    if abs(t(i+1)-round(t(i+1)/tscr)*tscr) < dt/2
        fprintf('\ntime = %7.2f (iter %5.1d)',t(i+1),i)
        
        figure(100); clf;
        subplot(4,2,1:4); surf(xx,zz,v(:,:,i),'EdgeColor','none');
                        colorbar('NO'); colormap(redblue)
                        cax = caxis; caxis([-1 1]*max(abs(cax)));
                        axis image; view(2);
                        xlabel('x'), ylabel('z'); title(sprintf('v(x,z,t = %.1f)',t(i+1)))
        subplot(4,2,5); plot(t(1:i+1),d(:,1:i+1));
                        ax = axis; axis([0 tend ax(3:4)]); grid on
                        xlabel('t'), ylabel('d(t)');
        subplot(4,2,7); plot(t(1:i+1),u(:,1:i+1));
                        ax = axis; axis([0 tend ax(3:4)]); grid on
                        xlabel('t'), ylabel('u(t)');
        subplot(4,2,6); plot(t(1:i+1),y(:,1:i+1)); 
                        ax = axis; axis([0 tend ax(3:4)]); grid on
                        xlabel('t'), ylabel('y(t)');
        subplot(4,2,8); plot(t(1:i+1),z(:,1:i+1)); 
                        ax = axis; axis([0 tend ax(3:4)]); grid on
                        xlabel('t'), ylabel('z(t)');
        drawnow
    end
                                         
    
end

fprintf(' END.\n')
