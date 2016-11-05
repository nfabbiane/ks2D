function [C,Cphys] = ks_init_output(pos,sigma,xx,zz)
%
%   [C,Cphys] = KS_init_input(pos,sigma,xx,zz)
%

nout = size(pos,1);
[nx,nz] = size(xx);
nq  = nx*nz;

% - compute input matrix in physical space
Cphys = zeros(nx,nz,nout);
for i = 1:nout
    Cphys(:,:,i)  = conj(fft2(exp(- ((xx-pos(i,1)).^2)/sigma(i,1).^2 ...
                                  - ((zz-pos(i,2)).^2)/sigma(i,2).^2 )/sqrt(prod(sigma(i,:)))));
end

% - reorder for state space formulation
C = zeros(nout,nq);
for i = 1:nout
    C(i,:) = reshape(Cphys(:,:,i),nq,1);
end