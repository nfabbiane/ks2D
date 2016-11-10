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
dt   = .25;         % time-step

% control parameters
rho     = 1e0;
maxiter = 1e1;
step    = 1e-5;




%% Initialization
[A,xx,zz,L,~,~,lambda] = ks_init(P,R,S,V,LX,LZ,Lf,NX,NZ);




%% Inputs matrix B

% disturbance d (Gaussian shape at x_d, z_d with sigma_d variance)
nd = 3; 
posd = zeros(nd,2); posd(:,1) = 0;
                    posd(:,2) = -LZ/2 + LZ/(2*nd):LZ/(nd):LZ/2 - LZ/(2*nd);
sigd = zeros(nd,2); sigd(:,1) = 4;
                    sigd(:,2) = 4;

[Bd,Bdfou] = ks_init_input(posd,sigd,xx,zz);

      
% actuator u (Gaussian shape at x_u, z_u with sigma_u variance)
nu = 3; 
posu = zeros(nu,2); posu(:,1) = 200;
                    posu(:,2) = -LZ/2 + LZ/(2*nu):LZ/(nu):LZ/2 - LZ/(2*nu);
sigu = zeros(nu,2); sigu(:,1) = 4;
                    sigu(:,2) = 4;

[Bu,Bufou] = ks_init_input(posu,sigu,xx,zz);

% control penalty (for control optimisation)
W = eye(nu) * rho;




%% Outputs matrix C

% measurement y (Gaussian shape at x_y, z_y with sigma_y variance)
ny = nu; 
posy = zeros(ny,2); posy(:,1) = 100;
                    posy(:,2) = -LZ/2 + LZ/(2*ny):LZ/(ny):LZ/2 - LZ/(2*ny);
sigy = zeros(ny,2); sigy(:,1) = 4;
                    sigy(:,2) = 4;

[Cy,Cyfou] = ks_init_output(posy,sigy,xx,zz,LX,LZ);


% output z (Gaussian shape at x_z with sigma_z variance)
nz = nu; 
posz = zeros(nz,2); posz(:,1) = 300;
                    posz(:,2) = -LZ/2 + LZ/(2*nz):LZ/(nz):LZ/2 - LZ/(2*nz);
sigz = zeros(nz,2); sigz(:,1) = 4;
                    sigz(:,2) = 4;

[Cz,Czfou] = ks_init_output(posz,sigz,xx,zz,LX,LZ);




%% Riccati-based optimisation
fprintf('\nKS Riccati-based optimal control. ')

% - solve Riccati eq.
X = care(A,Bu,Cz'*Cz,W);

% - compute control gains
KRic = -W\Bu'*X;

fprintf(' END.\n')




%% Adjoint-based optimisation
t = 0:dt:tend; nt = length(t);
nq = size(A,1);

% init optimisation variables
KAdj = zeros(NX,NZ,nu,maxiter);
JAdj = zeros(1,maxiter);
dJdK = zeros(NX,NZ,nu,maxiter);

% optimisation loop
fprintf('\nKS adjoint-based optimal control.\n')

for iter = 1:maxiter
    
    fprintf('\niter %2.0d: ',iter); tic
    
    % loop on disturbances
    for m = 1:nd

        % direct loop
        q = zeros(NX,NZ,nt); q(:,:,1) = Bdfou(:,:,m);
        f = zeros(NX,NZ,1);
        z = zeros(nz,nt);
        u = zeros(nu,nt);
        
        for i = 1:nt-1
            
            % output(s)
            for r = 1:nz
                z(r,i) = real(sum(sum(Czfou(:,:,r) .* q(:,:,i)))) * (LX*LZ);
            end
            
            % update cost function
            JAdj(iter) = JAdj(iter) + 1/2 * (z(:,i)' * z(:,i)) * dt;
            
            % control signal
            for r = 1:nu
                u(r,i) = real(sum(sum(KAdj(:,:,r,iter) .* q(:,:,i)))) * (LX*LZ);
            end
            
            % forcing
            f(:,:) = 0;
            for r = 1:nu
                f(:,:) = f(:,:) + Bufou(:,:,r) * u(r,i);
            end

            % KS time-step
            q(:,:,i+1) = ks_timestep(q(:,:,i),f,dt,L,lambda);
            
        end

        
        % adjoint loop
        l = zeros(NX,NZ,nt); l(:,:,end) = q(:,:,end);
        f = zeros(NX,NZ,1);
        w = zeros(nu,nt);
        
        for i = nt:-1:2
            
            % control signal
            for r = 1:nu
                w(r,i) = real(sum(sum(conj(Bufou(:,:,r)) .* l(:,:,i)))) * (LX*LZ);
            end
         
            % update gradient
            for r = 1:nu
                dJdK(:,:,r,iter) = dJdK(:,:,r,iter) + ...
                           (W(r,:)*u(:,i) - w(r,i)) .* conj(q(:,:,i)) * dt;
            end

            % forcing
            f(:,:) = 0;
            for r = 1:nz
                f(:,:) = f(:,:) - conj(Czfou(:,:,r)) * z(r,i);
            end
            for r = 1:nu
                f(:,:) = f(:,:) - conj(KAdj(:,:,r,iter)) * (W(r,:) * u(:,i));
            end
            for r = 1:nu
                f(:,:) = f(:,:) + conj(KAdj(:,:,r,iter)) * w(r,i);
            end

            % KS time-step
            l(:,:,i-1) = ks_timestep(l(:,:,i),f,dt,conj(L),lambda);

        end
    
    end
    
    
    % update control gains
    KAdj(:,:,:,iter+1) = KAdj(:,:,:,iter) - dJdK(:,:,:,iter) * step;
    
    
    % optimisation status
    runtime = toc;
    fprintf('J = %8.2e, |dJdK|_2 = %8.2e (runtime %.2fs)',...
             JAdj(iter),sum(sum(sum(dJdK(:,:,:,iter).*conj(dJdK(:,:,:,iter))))),runtime)

    figure(100); clf;
    subplot(5,1,1:2); surf(xx,zz,ifft2(KAdj(:,:,1,iter))*(NX*NZ),'EdgeColor','none');
                    colorbar('EO'); colormap(redblue)
                    cax = caxis; caxis([-1 1]*max(abs(cax)));
                    axis image; view(2);
                    xlabel('x'), ylabel('z'); title('K_1')
    subplot(5,1,3:4); surf(xx,zz,ifft2(dJdK(:,:,1,iter))*(NX*NZ),'EdgeColor','none');
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
    subplot(nu,2,1+2*(m-1)); surf(xx,zz,q2v(KRic(m,:).',NX,NZ)/(LX*LZ),'EdgeColor','none');
                    colorbar('EO'); colormap(redblue)
                    cax = caxis; caxis([-1 1]*max(abs(cax)));
                    axis image; view(2);
                    xlabel('x'), ylabel('z'); title(sprintf('K_%d Riccati-based',m))
    subplot(nu,2,2+2*(m-1)); surf(xx,zz,ifft2(KAdj(:,:,m,end))*(NX*NZ),'EdgeColor','none');
                    colorbar('EO'); colormap(redblue)
                    cax = caxis; caxis([-1 1]*max(abs(cax)));
                    axis image; view(2);
                    xlabel('x'), ylabel('z'); title(sprintf('K_%d adjoint-based',m))
end