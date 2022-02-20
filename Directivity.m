function [ Dir, Prad ] = Directivity( E_tot_abs, TH, dth, dph, er, R_FF )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%   Note: fix documentation of Directivity routine
    %% Constants
    e0 = 8.8541878128e-12;      % Permittivity of free space [F/m]
    u0 = 4*pi*1e-7;             % Permeability of free space [H/m]
    %% Calculate wave impedance
    Z = sqrt(u0 / (e0*er));
    %% Calculate radiated power
    Prad = 2 * sum( sum( E_tot_abs .* sin(TH) )) * dth * dph;
    %% Calculate radiation intensity
    U = (E_tot_abs.^2) .* (R_FF.^2) / (2 * Z);
    %% Calculate directivity
    Dir = U / (Prad / (4*pi));
end