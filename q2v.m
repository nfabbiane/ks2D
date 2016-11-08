function v = q2v(q,nx,nz)
%
%   v = q2v(q)
%

nqx = ceil(nx/2);
nqz = ceil(nz/2-1)*2+1;
nq  = nqx*nqz;

if nq == length(q)
    v = zeros(nx,nz);
    
    % figure(1000); subplot(1,3,1); imagesc(abs(reshape(q,nqx,nqz))); axis image; colorbar; size(v)
    
    ii = 1:nqx;
    jj = [1:(nqz+1)/2, (nqz+1)/2+1+(1:(nqz/2))];
    v(ii,jj) = reshape(q,nqx,nqz);
    
    % figure(1000); subplot(1,3,2); imagesc(abs(v)); axis image; colorbar; size(v)
    
    ii = nqx+(nqx:-1:2); vconj = conj(reshape(q,nqx,nqz));
    v(ii,jj) = vconj(2:end,:);
    
    % figure(1000); subplot(1,3,3); imagesc(abs(v)); axis image; colorbar; size(v)
    
    v = ifft2(v, 'symmetric') * nx*nz;
else
    v = [];
    disp('Error: nq must be equal to nx * nz!'); return
end