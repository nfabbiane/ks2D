clc, clear all, %close all
addpath('./lib')

%% Parameters

% Kuramoto
omegarmax = 2.5e-3; % maximum growth rate
alphamax = 1.5e-1;  % alpha of maximum growth rate
betamax = 1.5e-1;   % maximum unstable beta

V = 1.0;            % phase speed (= group speed, KS is not dispersive)

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
tend = 750;         % final time
dt   = 1.0;         % time-step

% control parameters
rho     = 1e0;
maxiter = 1e1;
step    = 1e0;




%% Initialization
[A,xx,zz] = ks_init(P,R,S,V,LX,LZ,Lf,NX,NZ);




%% Inputs matrix B

% disturbance d (Gaussian shape at x_d, z_d with sigma_d variance)
nd = 3; 
posd = zeros(nd,2); posd(:,1) = 0;
                    posd(:,2) = -LZ/2 + LZ/(2*nd):LZ/(nd):LZ/2 - LZ/(2*nd);
sigd = zeros(nd,2); sigd(:,1) = 4;
                    sigd(:,2) = 4;

Bd = ks_init_input(posd,sigd,xx,zz);

      
% actuator u (Gaussian shape at x_u, z_u with sigma_u variance)
nu = 3; 
posu = zeros(nu,2); posu(:,1) = 200;
                    posu(:,2) = -LZ/2 + LZ/(2*nu):LZ/(nu):LZ/2 - LZ/(2*nu);
sigu = zeros(nu,2); sigu(:,1) = 4;
                    sigu(:,2) = 4;

Bu = ks_init_input(posu,sigu,xx,zz);

% control penalty (for control optimisation)
W = eye(nu) * rho;




%% Outputs matrix C

% measurement y (Gaussian shape at x_y, z_y with sigma_y variance)
ny = nu; 
posy = zeros(ny,2); posy(:,1) = 100;
                    posy(:,2) = -LZ/2 + LZ/(2*ny):LZ/(ny):LZ/2 - LZ/(2*ny);
sigy = zeros(ny,2); sigy(:,1) = 4;
                    sigy(:,2) = 4;

Cy = ks_init_output(posy,sigy,xx,zz,LX,LZ);


% output z (Gaussian shape at x_z with sigma_z variance)
nz = nu; 
posz = zeros(nz,2); posz(:,1) = 300;
                    posz(:,2) = -LZ/2 + LZ/(2*nz):LZ/(nz):LZ/2 - LZ/(2*nz);
sigz = zeros(nz,2); sigz(:,1) = 4;
                    sigz(:,2) = 4;

Cz = ks_init_output(posz,sigz,xx,zz,LX,LZ);




%% Riccati-based optimisation
fprintf('\nKS Riccati-based optimisation. ')

% - solve Riccati eq.
X = care(A,Bu,Cz'*Cz,W);

% - compute control gains
KRic = -W\Bu'*X;

fprintf(' END.\n')




%% Adjoint-based optimisation
t = 0:dt:tend; nt = length(t);
nq = size(A,1);

% init optimisation variables
KAdj = zeros(nu,nq,maxiter);
JAdj = zeros(1,maxiter);
dJdK = zeros(nu,nq,maxiter);

% optimisation loop
fprintf('\nKS adjoint-based optimisation.\n')

for iter = 1:maxiter

    % controlled system
    Actr = sparse(A + Bu*KAdj(:,:,iter));
    Q = Cz'*Cz + KAdj(:,:,iter)'*W*KAdj(:,:,iter);
    
    % direct time-stepper
    Adtdir = sparse( (eye(size(Actr)) - Actr * dt/2) \ ...
                     (eye(size(Actr)) + Actr * dt/2) );
              
    % adjoint time-stepper
    Adtadj = sparse( (eye(size(Actr)) - Actr' * dt/2) \ ...
                     (eye(size(Actr)) + Actr' * dt/2) );
              
    
    % loop on disturbances
    for m = 1:nd
        
        % init state space variables
        q = zeros(nq,nt);
        l = zeros(nq,nt);
        f = zeros(nq,1);

        % direct loop
        q(:,1) = Bd(:,m);
        for i = 1:nt-1

            % forcing
            f(:) = 0;

            % KS time-step
            q(:,i+1) = Adtdir * (q(:,i) + f*dt);

            % compute cost function
            JAdj(iter) = JAdj(iter) + 1/2 * (q(:,i)'*Q*q(:,i)) * dt;

        end

        % adjoint loop
        l(:,end) = q(:,end);
        for i = nt:-1:2

            % forcing
            f(:) = -Q*q(:,i);

            % KS time-step
            l(:,i-1) = Adtadj * (l(:,i) + f*dt);

            % compute gradient
            dJdK(:,:,iter) = dJdK(:,:,iter) + ...
                     (W*KAdj(:,:,iter)*q(:,i) - Bu'*l(:,i)) * q(:,i)' * dt;

        end
    
    end
    
    
    % update control gains
    KAdj(:,:,iter+1) = KAdj(:,:,iter) - dJdK(:,:,iter) * step;
    
    
    % optimisation status
    fprintf('\niter %3.0d: J = %8.2e, |dJdK|_2 = %8.2e',...
                    iter,      JAdj(iter),       sum(sum(dJdK(:,:,iter)*dJdK(:,:,iter)')))

    figure(100); clf;
    subplot(5,1,1:2); surf(xx,zz,q2v(KAdj(1,:,iter).',NX,NZ),'EdgeColor','none');
                    colorbar('EO'); colormap(redblue)
                    cax = caxis; caxis([-1 1]*max(abs(cax)));
                    axis image; view(2);
                    xlabel('x'), ylabel('z'); title('K_1')
    subplot(5,1,3:4); surf(xx,zz,q2v(dJdK(1,:,iter).',NX,NZ),'EdgeColor','none');
                    colorbar('EO'); colormap(redblue)
                    cax = caxis; caxis([-1 1]*max(abs(cax)));
                    axis image; view(2);
                    xlabel('x'), ylabel('z'); title('dJ/dK_1')
    subplot(5,1,5); semilogy(1:maxiter,JAdj,'s-');
                    ax = axis; axis([1 maxiter ax(3:4)]); grid on
                    xlabel('iter'), ylabel('J');
	drawnow
    
end

fprintf(' END.\n')




%% Compare
figure(1); clf
for m = 1:nu
    subplot(nu,2,1+2*(m-1)); surf(xx,zz,q2v(KRic(m,:).',NX,NZ),'EdgeColor','none');
                    colorbar('EO'); colormap(redblue)
                    cax = caxis; caxis([-1 1]*max(abs(cax)));
                    axis image; view(2);
                    xlabel('x'), ylabel('z'); title(sprintf('K_%d Riccati',m))
    subplot(nu,2,2+2*(m-1)); surf(xx,zz,q2v(KAdj(m,:,end).',NX,NZ),'EdgeColor','none');
                    colorbar('EO'); colormap(redblue)
                    cax = caxis; caxis([-1 1]*max(abs(cax)));
                    axis image; view(2);
                    xlabel('x'), ylabel('z'); title(sprintf('K_%d Adjoint',m))
end