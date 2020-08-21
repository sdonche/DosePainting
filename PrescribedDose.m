function PresDose = PrescribedDose(I,I_high,I_low,D_high,D_low)
% PrescribedDose calculates the prescribed dose from a (PET) intensity.
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
% D(high)       = 28 Gy
% D(low)        = 20 Gy
% I(high)       = 95% of PET voxel intensity        
% I(low)        = I(high)*0.25

PresDose = D_low + (I - I_low)/(I_high - I_low)*(D_high - D_low);

end

