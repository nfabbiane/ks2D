function [B,Bfou,Bphy] = ks_init_input(pos,sigma,xx,zz)
%
%   [B,Bphys,Bphy] = KS_init_input(pos,sigma,xx,zz)
%

nin = size(pos,1);
[nx,nz] = size(xx);
nq  = nx*nz;

% - compute input matrix in physical space
Bphy = zeros(nx,nz,nin);
Bfou = zeros(nx,nz,nin);

for l = 1:nin
    Bphy(:,:,l)  = exp(- ((xx-pos(l,1)).^2)/sigma(l,1).^2 ...
                       - ((zz-pos(l,2)).^2)/sigma(l,2).^2 )/sqrt(prod(sigma(l,:)));
                                   
    Bfou(:,:,l)  = fft2(Bphy(:,:,l)) / (nx*nz);
end

% - reorder for state space formulation
B = zeros(nq,nin);

for l = 1:nin
    B(:,l) = v2q(Bphy(:,:,l));
end