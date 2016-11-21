function [B,Bfou,Bphy,dB,dBfou,dBphy] = ks_init_input(pos,sigma,xx,zz,LX,LZ)
%
%   [B,Bphys,Bphy] = KS_init_input(pos,sigma,xx,zz)
%

nin = size(pos,1);
[nx,nz] = size(xx);
nq  = nx*nz;

% - compute input matrix in physical space
Bphy  = zeros(nx,nz,nin);
Bfou  = zeros(nx,nz,nin);
dBphy = zeros(nx,nz,nin,2);
dBfou = zeros(nx,nz,nin,2);

for l = 1:nin
    arg = - ((mod(xx-pos(l,1)     ,LX)     ).^2)/sigma(l,1).^2 ...
          - ((mod(zz-pos(l,2)+LZ/2,LZ)-LZ/2).^2)/sigma(l,2).^2 ;
      
    Bphy(:,:,l) = exp(arg)/sqrt(prod(sigma(l,:)));
                   
    dBphy(:,:,l,1)  = -1/sigma(l,1) ...
                      * (-2 * (mod(xx-pos(l,1)     ,LX)     )/sigma(l,1)) ...
                     .* arg .* Bphy(:,:,l);
    dBphy(:,:,l,2)  = -1/sigma(l,2) ...
                      * (-2 * (mod(zz-pos(l,2)+LZ/2,LZ)-LZ/2)/sigma(l,2)) ...
                     .* arg .* Bphy(:,:,l);
                                   
    Bfou(:,:,l)    = fft2( Bphy(:,:,l)  ) / (nx*nz);
    dBfou(:,:,l,1) = fft2(dBphy(:,:,l,1)) / (nx*nz);
    dBfou(:,:,l,2) = fft2(dBphy(:,:,l,2)) / (nx*nz);
end

% - reorder for state space formulation
B  = zeros(nq,nin);
dB = zeros(nq,nin,2);

for l = 1:nin
    B(:,l)    = v2q( Bphy(:,:,l)  );
    dB(:,l,1) = v2q(dBphy(:,:,l,1));
    dB(:,l,2) = v2q(dBphy(:,:,l,2));
end