function [B,Bphys] = ks_init_input(pos,sigma,xx,zz)
%
%   [B,Bphys] = KS_init_input(pos,sigma,xx,zz)
%

nin = size(pos,1);
[nx,nz] = size(xx);
nq  = nx*nz;

% - compute input matrix in physical space
Bphys = zeros(nx,nz,nin);
for i = 1:nin
    Bphys(:,:,i)  = 1/(nx*nz)*fft2(exp(- ((xx-pos(i,1)).^2)/sigma(i,1).^2 ...
                                       - ((zz-pos(i,2)).^2)/sigma(i,2).^2 )/sqrt(prod(sigma(i,:))));
end

% - reorder for state space formulation
B = zeros(nq,nin);
for i = 1:nin
    B(:,i) = reshape(Bphys(:,:,i),nq,1);
end