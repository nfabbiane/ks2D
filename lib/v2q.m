function q = v2q(v)
%
%   q = v2q(v)
%

[nx,nz] = size(v);
nq = nx*nz;

q = reshape(fft2(v)/nq,nq,1);