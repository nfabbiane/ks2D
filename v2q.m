function q = v2q(v)
%
%   q = v2q(v)
%

[nx,nz] = size(v);

nqx = ceil(nx/2);
nqz = ceil(nz/2-1)*2+1;
nq  = nqx*nqz;

ii = 1:nqx;
jj = [1:(nqz+1)/2, (nqz+1)/2+1+(1:(nqz/2))];

v = fft2(v) / (nx*nz);
q = reshape(v(ii,jj),nq,1);