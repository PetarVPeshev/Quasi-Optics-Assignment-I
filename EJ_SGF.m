function [ ej_SGF ] = EJ_SGF( er, k, kx, ky )
%EJ_SGF This function computes the full Spectral Green's Function (SGF)
%   The length of the x and y-component of the wave number are extracted
%   from the arrays to be used later in partitioning the SGF matrix. In
%   addition, based on the provided relative permittivity, the wave
%   impedance is calculated to be used in the SGF derivation. Moreover, the
%   z-component wave number is calculated from the provided x and
%   y-components and the magnitude of the wave number. Then, for the SGF, a
%   matrix is computed based on the cartesian coordinate wave number
%   components, and then partitioned and multiplied by the respective
%   constant for each x and y-component of the wave number. The output
%   matrix is 4D with dimension one representing the different x-component
%   values of the wave number, dimension two the y-component values,
%   dimension three the rows of the SGF matrix, and dimension four the
%   columns of the SGF matrix.
%   Note: fix documentation of EJ_SGF routine
    %% Extract wave number x and y-component length
    kx_l = length(kx);
    ky_l = length(ky);
    %% Create meshgrid
    [ KX, KY ] = meshgrid(kx, ky);
    %% Constants
    e0 = 8.8541878128e-12;      % Permittivity of free space [F/m]
    u0 = 4*pi*1e-7;             % Permeability of free space [H/m]
    %% Calculate wave impedance
    Z = sqrt(u0 / (e0*er));
    %% Calculate z-component of wave number
    kz = -1j*sqrt(-(k.^2 - KX.^2 - KY.^2));
    %% Calcualte Spectral Green's Functions (SGF)
    G = [(k^2 - KX.^2) (-KX .* KY) (-KX .* kz);
         (-KY .* KX) (k^2 - KY.^2) (-KY .* kz);
         (-kz .* KX) (-kz .* KY) (k^2 - kz.^2)];
    ej_SGF = -Z ./ (2 * k .* kz) .* ones(kx_l, ky_l, 3, 3);
    for i = 1:3
        for n = 1:3
            ej_SGF(:, :, i, n) = ej_SGF(:, :, i, n) .* G(((i-1)*kx_l + 1) : (i*kx_l), ...
                                                         ((n-1)*ky_l + 1) : (n*ky_l));
        end
    end
end