function [ ej_SGF ] = EJ_SGF( er, k, kx, ky )
%EJ_SGF This function computes the full Spectral Green's Function (SGF)
%   The function takes the relative permittivity of the medium, the
%   magnitude of the wave number, and meshgrid of the x and y-components of
%   the wave number as inputs.
%   It outputs the Spectral Green's Function (SGF) of the electric field
%   due to electric current in a 4D matrix with dimensions 1 and 2
%   corresponding to the  x and y-component meshgrid of the wave number,
%   and dimensions 3 and 4 being the dyadic tensor product.
    %% Extract wave number x and y-component length
    kx_l = size(kx, 2);
    ky_l = size(ky, 1);
    %% Constants
    e0 = 8.8541878128e-12;      % Permittivity of free space [F/m]
    u0 = 4*pi*1e-7;             % Permeability of free space [H/m]
    %% Calculate wave impedance
    Z = sqrt(u0 / (e0*er));
    %% Calculate z-component of wave number
    kz = -1j * sqrt( -(k.^2 - kx.^2 - ky.^2) );
    %% Calcualte Spectral Green's Functions (SGF)
    G = [(k^2 - kx.^2) (-kx .* ky) (-kx .* kz);
         (-ky .* kx) (k^2 - ky.^2) (-ky .* kz);
         (-kz .* kx) (-kz .* ky) (k^2 - kz.^2)];
    ej_SGF = -Z ./ (2 * k .* kz) .* ones(kx_l, ky_l, 3, 3);
    for i = 1:3
        for n = 1:3
            ej_SGF(:, :, i, n) = ej_SGF(:, :, i, n) .* G(((i-1)*kx_l + 1) : (i*kx_l), ...
                                                         ((n-1)*ky_l + 1) : (n*ky_l));
        end
    end
end