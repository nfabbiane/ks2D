function [q,v] = ks_timestep(q,f,dt,A,lambda)

[nx,nz] = size(q);
nq  = nx*nz;

% Runge-Kutta coefficients
ns = 3;
a  = [0     1/3     3/4     1 ];
b  = [0    -5/9  -153/128     ];
c  = [1/3  15/16    8/15      ];

h = 0.*q;

% Runge-Kutta
for m = 1:ns
    dts=  dt*(a(m+1) - a(m));
    
    h = f + b(m)*h;
    q = ((1 + dts/2 * A) .* q + c(m)*dt * h)./(1 - dts/2 * A);
end

% Fringe forcing
v = (1 - lambda) .* (ifft2(q, 'symmetric') * nq);
q = fft2(v) / nq;