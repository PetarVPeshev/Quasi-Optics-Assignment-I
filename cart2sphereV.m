function [ Asph ] = cart2sphereV( A, TH, PH )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   Note: fix documentation of cart2sphere routine
    %% Extract length of theta and phi
    TH_l = size(TH, 2);
    PH_l = size(PH, 1);
    %% Define transformation matrix
    TMatrix = [(sin(TH) .* cos(PH))   (sin(TH) .* sin(PH))      cos(TH);
               (cos(TH) .* cos(PH))   (cos(TH) .* sin(PH))     -sin(TH);
                     -sin(PH)                cos(PH)        zeros(TH_l, PH_l)];
    %% Partition transformation matrix
    TM = zeros(TH_l, PH_l, 3, 3);
    for i = 1:3
        for n = 1:3
            TM(:, :, i, n) = TMatrix(((i-1)*TH_l + 1) : (i*TH_l), ...
                                     ((n-1)*PH_l + 1) : (n*PH_l));
        end
    end
    %% Reshape matricies
    TM = permute(TM, [3 4 1 2]);
    A = permute(A, [3 1 2]);
    A = reshape(A, [3 1 TH_l PH_l]);
    %% Convert to spherical coordinates
    Asph = pagemtimes(TM, A);
    %% Reshape output to standard dimensions
    Asph = permute(Asph, [3 4 1 2]);
end