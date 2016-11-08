function [C,Cfou,Cphy] = ks_init_output(pos,sigma,xx,zz)
%
%   [C,Cfou,Cphy] = KS_init_input(pos,sigma,xx,zz)
%

nout = size(pos,1);
[nx,nz] = size(xx);
nq  = nx*nz;

% - compute input matrix in physical space
Cphy = zeros(nx,nz,nout);
Cfou = zeros(nx,nz,nout);

for l = 1:nout
    Cphy(:,:,l)  = exp(- ((xx-pos(l,1)).^2)/sigma(l,1).^2 ...
                       - ((zz-pos(l,2)).^2)/sigma(l,2).^2 )/sqrt(prod(sigma(l,:)));
                   
    Cfou(:,:,l)  = conj(fft2(Cphy(:,:,l)));
end

% - reorder for state space formulation
C = zeros(nout,nq);

for l = 1:nout
    C(l,:) = v2q(Cphy(:,:,l)) * (nx*nz);
end