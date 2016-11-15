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
tend = 500;         % final time
dt   = 1.0;         % time-step

% control parameters
rho     = 1e4;      % control penalty
maxiter = 40;        % max number of iterations
epsilon = 1e-5;     % stop tollerance: |(J_i - J_i-1)/J_i-1| < eps




%% Initialization
[A,xx,zz,L,~,~,lambda] = ks_init(P,R,S,V,LX,LZ,Lf,NX,NZ);




%% Inputs matrix B

% disturbance d (Gaussian shape at x_d, z_d with sigma_d variance)
nd = 12; 
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

% initialize control gains
K = zeros(NX,NZ,nu,maxiter);
% for r = 1:nu
%     K(:,:,r,1) = fft2(q2v(KRic(r,:).',NX,NZ)/(LX*LZ))/(NX*NZ);
% end


% initialize gradient
dJdK  = zeros(NX,NZ,nu,maxiter);
dJdKn = zeros(1,maxiter);

% initialize cost function
J  = zeros(1,maxiter);
dJ = zeros(1,maxiter);

% initialize conjuagate gradient coefficients
P     = zeros(NX,NZ,nu);
dJdP  = zeros(NX,NZ,nu);
up    = zeros(nu,1);
alpha = zeros(1,maxiter);
beta  = zeros(1,maxiter);

% optimisation loop
fprintf('\nKS adjoint-based optimal control.\n')

for iter = 1:maxiter
    
    fprintf('\niter %2.0d: ',iter); tic
    
    % compute gradient and cost function
    q = zeros(NX,NZ,nt,nd);
    z = zeros(nz,nt,nd);
    u = zeros(nu,nt,nd);

    l = zeros(NX,NZ,nt,nd);
    h = zeros(nu,nt,nd);

    f = zeros(NX,NZ,1);
    
    % loop on disturbances
    for m = 1:nd

        % direct loop
        q(:,:,1,m) = Bdfou(:,:,m);
        
        for i = 1:nt-1
            
            % output(s)
            for r = 1:nz
                z(r,i,m) = real(sum(sum(Czfou(:,:,r) .* q(:,:,i,m)))) * (LX*LZ);
            end
            
            % control signal
            for r = 1:nu
                u(r,i,m) = real(sum(sum(K(:,:,r,iter) .* q(:,:,i,m)))) * (LX*LZ);
            end
            
            % update cost function
            J(iter) = J(iter) + 1/2 * (z(:,i,m)' * z(:,i,m) ...
                                          + u(:,i,m)' * W * u(:,i,m)) * dt;
            
            % forcing
            f(:,:) = 0;
            for r = 1:nu
                f(:,:) = f(:,:) + Bufou(:,:,r) * u(r,i,m);
            end

            % KS time-step
            q(:,:,i+1,m) = ks_timestep(q(:,:,i,m),f,dt,L,lambda);
            
        end

        % adjoint loop
        l(:,:,end,m) = q(:,:,end,m);
        
        for i = nt:-1:2
            
            % control signal
            for r = 1:nu
                h(r,i,m) = real(sum(sum(conj(Bufou(:,:,r)) .* l(:,:,i,m)))) * (LX*LZ);
            end
         
            % update gradient
            for r = 1:nu
                dJdK(:,:,r,iter) = dJdK(:,:,r,iter) + ...
                           (W(r,:)*u(:,i,m) - 1/2 * h(r,i,m)) .* conj(q(:,:,i,m)) * dt;
            end

            % forcing
            f(:,:) = 0;
            for r = 1:nz
                f(:,:) = f(:,:) - 2 * conj(Czfou(:,:,r)) * z(r,i,m);
            end
            for r = 1:nu
                f(:,:) = f(:,:) - 2 * conj(K(:,:,r,iter)) * (W(r,:) * u(:,i,m)) ...
                                +     conj(K(:,:,r,iter)) *  h(r,i,m);
            end

            % KS time-step
            l(:,:,i-1,m) = ks_timestep(l(:,:,i,m),f,dt,conj(L),lambda);

        end
    
    end
    
    
    % cost function variation
    if iter > 1
        dJ(iter) = abs(J(iter)-J(iter-1))/abs(J(iter-1));
    end
    
    
    % conjugate gradient
    
    % - direction
    dJdKn(iter) = sum(sum(sum(conj(dJdK(:,:,:,iter)).*dJdK(:,:,:,iter))));
    if iter > 1
        beta(iter) = dJdKn(iter) / dJdKn(iter-1); % Fletcher?Reeves
    end
    P(:,:,:) = - dJdK(:,:,:,iter) + beta(iter) * P;
    
    % - step length (solution of the approximated quadratic problem)
    dJdP(:,:,:) = 0;
    for m = 1:nd
        for i = 1:nt-1
            for r = 1:nu
                up(r) = real(sum(sum(P(:,:,r) .* q(:,:,i,m)))) * (LX*LZ);
            end
            for r = 1:nu
                dJdP(:,:,r) = dJdP(:,:,r) + W(r,:)*up .* conj(q(:,:,i,m)) * dt;
            end
        end
    end
    alpha(iter) = - sum(sum(sum(conj(dJdP) .* dJdK(:,:,:,iter)))) / ...
                                        sum(sum(sum(conj(dJdP) .* dJdP)));
    alpha(iter) = real(alpha(iter)); % numerical error: imag(alpha) ~ 1e-18
    
    % - update
    K(:,:,:,iter+1) = K(:,:,:,iter) + alpha(iter) * P;
    
    
    % visualize status
    runtime = toc;
    fprintf('J = %8.2e, |dJdK|_2 = %8.2e (runtime %.2fs)',...
                 J(iter),       dJdKn(iter),    runtime)

    figure(100); clf; iu = floor(nu/2)+1;
    subplot(6,1,1:2); surf(xx,zz,ifft2(K(:,:,iu,iter))*(NX*NZ),'EdgeColor','none');
                    hc = colorbar('EO'); colormap(redblue)
                    cax = caxis; caxis([-1 1]*max(abs(cax)));
                    axis image; view(2); shading interp
                    xlabel('x'), ylabel('z'); ylabel(hc,sprintf('K_%d',iu))
    subplot(6,1,3:4); surf(xx,zz,ifft2(dJdK(:,:,iu,iter))*(NX*NZ),'EdgeColor','none');
                    hc = colorbar('EO'); colormap(redblue)
                    cax = caxis; caxis([-1 1]*max(abs(cax)));
                    axis image; view(2); shading interp
                    xlabel('x'), ylabel('z'); ylabel(hc,sprintf('dJ/dK_%d',iu))
    subplot(6,1,5); semilogy(1:maxiter,dJ,'s-',[1 maxiter],[1 1]*epsilon,'-r');
                    ax = axis; axis([1 maxiter ax(3:4)]); grid on
                    xlabel('iter'), ylabel('\deltaJ');
    subplot(6,1,6); semilogy(1:maxiter,dJdKn,'s-');
                    ax = axis; axis([1 maxiter ax(3:4)]); grid on
                    xlabel('iter'), ylabel('||\partialJ/\partialK||_2');
	drawnow
    
    
    % stop condition
    if (iter > 1) & (dJ(iter) < epsilon); break; end
    
end

fprintf(' END.\n')




%% Compare
figure(1); clf
for m = 1:nu
    subplot(nu,2,1+2*(m-1)); surf(xx,zz,q2v(KRic(m,:).',NX,NZ)/(LX*LZ),'EdgeColor','none');
                    colorbar('EO'); colormap(redblue)
                    cax = caxis; caxis([-1 1]*max(abs(cax)));
                    axis image; view(2); shading interp
                    xlabel('x'), ylabel('z'); title(sprintf('K_%d Riccati-based',m))
    subplot(nu,2,2+2*(m-1)); surf(xx,zz,ifft2(K(:,:,m,iter)*(NX*NZ)),'EdgeColor','none');
                    colorbar('EO'); colormap(redblue)
                    cax = caxis; caxis([-1 1]*max(abs(cax)));
                    axis image; view(2); shading interp
                    xlabel('x'), ylabel('z'); title(sprintf('K_%d adjoint-based',m))
end