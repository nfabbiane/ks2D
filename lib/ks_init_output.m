function [C,Cfou,Cphy,dC,dCfou,dCphy] = ks_init_output(pos,sigma,xx,zz,LX,LZ)
%
%   [C,Cfou,Cphy] = KS_init_input(pos,sigma,xx,zz)
%

nout = size(pos,1);
[nx,nz] = size(xx);
nq  = nx*nz;

% - compute input matrix in physical space
Cphy = zeros(nx,nz,nout);
Cfou = zeros(nx,nz,nout);
dCphy = zeros(nx,nz,nout,2);
dCfou = zeros(nx,nz,nout,2);

for l = 1:nout
    arg = - (mod(xx-pos(l,1),LX).^2)/sigma(l,1).^2 ...
          - (mod(zz-pos(l,2),LZ).^2)/sigma(l,2).^2 ;
      
    Cphy(:,:,l) = exp(arg)/sqrt(prod(sigma(l,:)));
                   
    dCphy(:,:,l,1)  = -1/sigma(l,1) ...
                      * (-2 * (mod(xx-pos(l,1)     ,LX)     )/sigma(l,1)) ...
                     .* arg .* Cphy(:,:,l);
    dCphy(:,:,l,2)  = -1/sigma(l,2) ...
                      * (-2 * (mod(zz-pos(l,2)+LZ/2,LZ)-LZ/2)/sigma(l,2)) ...
                     .* arg .* Cphy(:,:,l);
                   
    Cfou(:,:,l)    = conj(fft2( Cphy(:,:,l)  )) / (nx*nz);
    dCfou(:,:,l,1) = conj(fft2(dCphy(:,:,l,1))) / (nx*nz);
    dCfou(:,:,l,2) = conj(fft2(dCphy(:,:,l,2))) / (nx*nz);
end

% - reorder for state space formulation
C  = zeros(nout,nq);
dC = zeros(nout,nq,2);

for l = 1:nout
    C(l,:)    = conj(v2q( Cphy(:,:,l)  )) * LX*LZ;
    dC(l,:,1) = conj(v2q(dCphy(:,:,l,1))) * LX*LZ;
    dC(l,:,2) = conj(v2q(dCphy(:,:,l,2))) * LX*LZ;
end