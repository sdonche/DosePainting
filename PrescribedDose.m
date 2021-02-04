function PresDose = PrescribedDose(I,I_high,I_low,D_high,D_low)
% PrescribedDose calculates the prescribed radiation dose pixelwise in function of image intensity.
% In this method, the dose increases linearly with intensity between D(low) 
% and D(high), respectively corresponding to I(low) and I(high).
% This method was adopted from [18F]FDG PET voxel intensity-based IMRT for
% head and neck cancer at Ghent University Hospital
% (https://doi.org/10.1016/j.radonc.2006.03.003)
% 
% Formula:
% D(I) = D(low) + (I - I(low))/(I(high) - I(low))*(D(high) - D(low))
% 
% with
% I             = voxel intensity
% D(I)          = prescribed dose
% D(high)       = maximal dose (e.g. 28 Gy)
% D(low)        = minimal dose (e.g. 20 Gy)
% I(high)       = 95% of PET voxel intensity        
% I(low)        = I(high)*0.25
%
% Written by Sam Donche

    if I >= I_high              % Cap off the prescribed dose at D(high)
        PresDose = D_high;
        
    elseif I <= I_low           % Minimal dose of D(low) 
        PresDose = D_low;
        
    else                        % Linear function between I(low) and I(high)
        PresDose = D_low + (I - I_low)/(I_high - I_low)*(D_high - D_low);
    
    end
end

