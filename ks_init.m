function [A,xx,zz,L,alpha,beta,lambda,LAMBDA] = ks_init(P,R,S,V,Lx,Lz,Lf,nx,nz)
%
% [A,x,z,L,alpha,beta,lambda,LAMBDA] = ks_init(P,R,S,V,Lx,Lz,Lf,nqx,nqz)
%
%   Linearised Kuramoto-Sivashinsky equation for v(x,z,t).
% 


% Physical grid
dx = Lx / (nx);
x  = (0:nx-1)' * dx;

dz = Lz / (nz);
z  = (0:nz-1)' * dz - Lz/2;

[zz,xx] = meshgrid(z,x);


% Linear operator
% - wavenumbers
alpha0 = 2*pi / Lx;
ialpha = ifftshift(-floor(nx/2):ceil(nx/2-1))';
alpha  = alpha0 * ialpha;

beta0 = 2*pi / Lz;
ibeta = ifftshift(-floor(nz/2):ceil(nz/2-1))';
beta  = beta0 * ibeta;

[beta,alpha] = meshgrid(beta,alpha);

% - modified KS equation
L = -V*1i*alpha.*(1 + beta.^2/(8*P)) + (P*alpha.^2 - (alpha.^4 + S*beta.^4))/R;
%L = -V*1i*alpha + (P*(alpha.^2+beta.^2) - (alpha.^4 + beta.^4))/R;
%L = -V*1i*alpha + (P*alpha.^2 - (alpha.^4 + beta.^4))/R;
%L = -V*1i*alpha.*(1 + beta.^2/(8*P)) + (P*alpha.^2 - (alpha.^4 + beta.^4))/R;

% - fringe (physical space)
Ls    = Lx - Lf;
lrise = .60 * Lf;
lfall = .40 * Lf;

lambda  = 0.8 * ( F((xx-Ls)/lrise) - F((xx-Lx)/lfall + 1) );


% Linear operator (state-space system)
% - reorder modes in a vector
nq = nx*nz;
l = reshape(L,nq,1);
A = spdiags(l,0,nq,nq);

% - fringe
lambdaf = fft(lambda(:,1))/nx;
LAMBDA  = sparse(zeros(nq,nq));

LAMBDAmod = zeros(nx);
for i = ialpha'
    for m = ialpha'
        LAMBDAmod((ialpha' == i),(ialpha' == i-m)) = lambdaf(ialpha == m);
    end
end

for j = 1:nz
    ii = (1:nx) + (j-1)*nx;
    LAMBDA(ii,ii) = LAMBDAmod;
end

% - full operator
A = sparse(A-LAMBDA);

end




%%% SUPPORT FUNCTION (fringe)
function s = F(xx)
    s = 0.*xx;
    sel = 0<xx & xx<1; s(sel) = 1./(1+exp(1./(xx(sel)-1)+1./xx(sel)));
    sel = xx>=1;       s(sel) = 1;
end