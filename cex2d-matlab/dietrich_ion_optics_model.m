close all; clear all; clc

%% Constants 
m = 131*(1e-3)/(6.022e23);
q = 1.66e-19;
eps0 = 8.854e-12;


%% Main
L = 0.5;
z     = linspace(0, L);
alpha = 0.5;
phi0  = 1;
phi_k = 5;
N     = 5;

G     = phi_k / phi0;
g     = G.^(1/N);
phi   = phi0 + phi_k * z.^(1/alpha);
zi = L;
r1 = 0.3;
P = I / (phi0 - phi_k).^(3/2)
P0 = 2 * P;
thetap = 0.625 * (r1/z1) * ((P/P0) - 1) * (1 - 0.195 * (P/P0 - 1))
for i=1:N

    dzi = Dzi(L, g, G, alpha, i)
end

%plot(z, phi)



%% Functions 
function dzi = Dzi(L, g, G, alpha, i)
% Dietrich
    dzi = L * ( (g^alpha - 1) / (G^alpha - 1) ) * g^(alpha*(i-1));
end

function k = calc_K(g)
    k = 2 * (1 - (sqrt(g)-1)/(g-1));
end

function Pi = perveance_i(I, phi_i, phi_im1)
    Pi = I / ( (phi_im1 - phi_i) / 2 )^(3/2);
end


function ri_zi = segment_radius(r_im1, dr_im1, ni, dzi, k, Pi, m, q, eps0)
    term2 = dr_im1*ni*dzi*k;
    term3 = (dzi.^2 / (r_im1)) * sqrt(m/(2*q)) * Pi/(2 * pi * (2.09.^2)*eps0);
    ri_zi = r_im1 + term2 + term3;

end




