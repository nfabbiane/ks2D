function v = q2v(q,nx,nz)
%
%   v = q2v(q)
%

nq = length(q);

if nq == nx*nz
    v = ifft2(reshape(q,nx,nz) * nq, 'symmetric');
else
    v = [];
    disp('Error: nq must be equal to nx * nz!'); return
end