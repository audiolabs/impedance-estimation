close all; clear; clc

% An example of the impedance estimator for a simulated impedance tube is 
% provided. It is assumed that the eigenvalues have been accurately 
% estimated. A finite element model with linear elements is used to 
% estimate the impedance of two samples. The a priori eigenvectors 
% are given by the analytical description of modes in a 1D rigid-walled 
% duct.

c0 = 344;           % speed of sound
L = 1.038;          % tube length
zref = [2-2i 2+2i]; % reference impedances
%---------------------------%
% estimated eigenvalues
evals = 1e2.*[1.5146+0.1303i 3.1720+0.1304i 4.8302+0.1305i;
    1.7996+0.1303i 3.4572+0.1304i 5.1156+0.1305i];
%---------------------------%
% impedance estimates
for zz = 1:length(zref)
    N = 100;    % number of finite elements
    h = L/N;    % element size
    x = 0:h:L;  % 1D mesh
    K = zeros(N+1, N+1);
    M = zeros(N+1, N+1);
    C = zeros(N+1, N+1);
    Ke = [1 -1 ; -1 1];
    Me = [1/3 1/6; 1/6 1/3];
    for ii = 1:N
        K(ii:ii+1, ii:ii+1) = K(ii:ii+1, ii:ii+1) + Ke;
        M(ii:ii+1, ii:ii+1) = M(ii:ii+1, ii:ii+1) + Me;
    end
    K = 1/h*K;          % stiffness matrix
    M = h/c0^2*M;       % mass matrix
    C(end, end) = L/c0; % damping matrix
    
    vec0 = zeros(length(x), 3); % eigenvectors of system (0)
    for jj = 1:3
        vec0(:, jj) = cos (jj*pi*x./L);
    end
    %
    zeta = get_imp_est(evals(zz, :), vec0, K, M, C);
    error = 100*abs(zref(zz)-zeta)./abs(zref(zz));
    disp(['impedance estimation error, in %: ' num2str(error)])
end
%---------------------------%
function zeta = get_imp_est(evals, vec0, K, M, C1)
Q = 5;  % number of iterations
for mm = 1:Q
    if mm == 1
        zeta = zeros(size(evals)) + 1e6;
    else
        zeta = zest;
    end
    zest = evals*0;
    for kk = 1:length(evals)
        C0 = 1/zeta(kk)*C1;
        [vec, val] = get_vec_val(vec0(:, kk), evals(kk), K, M, C0);
        lam0 = 2*pi*val;
        lam1 = 2*pi*evals(kk);
        S0 = lam0*vec.'*(C0 + 2i*lam0*M)*vec;
        S1 = lam1*vec.'*(2i*lam0*M)*vec;
        S2 = lam1*vec.'*C1*vec;  % Eq. (12)
        %
        zest(kk) = S2./(S0-S1); % Eq. (11)
    end
end
end
%---------------------------%
function [vec, val] = get_vec_val(vec, val, K, M, C)
%
lam = 2i*pi*val;
S = sparse(K + lam*C + lam^2*M);
sol = S\vec;
vec =  sol./norm(sol); % Eq. (13)
lam1 = ( -vec.'*C*vec + sqrt((vec.'*C*vec)^2 - ...
    4*(vec.'*M*vec)*(vec.'*K*vec)) )/(2*vec.'*M*vec);
lam2 = ( -vec.'*C*vec - sqrt((vec.'*C*vec)^2 - ...
    4*(vec.'*M*vec)*(vec.'*K*vec)) )/(2*vec.'*M*vec);
val1 = lam1/(2i*pi);
val2 = lam2/(2i*pi);
val1 = abs(real(val1)) + 1i*abs(imag(val1));
val2 = abs(real(val2)) + 1i*abs(imag(val2));
vals = [val1 val2];
tmp = min(abs(real(val)-real(vals)))== ...
    abs(real(val)-real(vals));
val = vals(tmp);
end

% publish('Appendix.m', struct('format','latex','outputDir','ltx-src'));