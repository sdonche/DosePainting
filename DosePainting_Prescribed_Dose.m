%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prescribed Dose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Information
% 
% D(I) = D(low) + (I - I(low))/(I(high) - I(low))*(D(high) - D(low))
% 
% with
% I             = voxel intensity
% D(I)          = prescribed dose
% D(high)       = 28 Gy
% D(low)        = 20 Gy
% I(high)       = 95% of PET voxel intensity        
% I(low)        = I(high)*0.25
%
% To minimize the influence of PET signal noise on I(low) and I(high), the
% dose was escalated between 25 and 100% of the 95th percentile PET voxel
% intensity value within PTV(69+PET)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLEAN SLATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;
clc;
imtool close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB TOOLBOXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/download/)
%   NIFTI and ANALYZE tools (https://nl.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
%   Medical Image Reader and Viewer (https://www.mathworks.com/matlabcentral/fileexchange/53745-medical-image-reader-and-viewer)
%   nifti_utils (https://github.com/justinblaber/nifti_utils)

% TODO update this in all scripts
disp('Reading Toolboxes...')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Extra\spm12')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Extra\NIfTI_20140122')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Extra\Medical Image Reader and Viewer')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Extra\nifti_utils-master\nifti_utils') % "updated" version of NIfTI_20140122
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Preclinical\Scripts\App_BoundingBox')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Preclinical\DosePainting\')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ IMAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Workdirectory
cd 'D:\Desktop\PhD Meeting\Demo\';
pathname = [pwd,'\'];

% Coregistered Images
pathname_coreg = [pathname, '5Coregister\'];
scans_coreg = dir(pathname_coreg);
scans_coreg = scans_coreg(3:end);

% Load coregistered images
disp(['Reading coregistered images from ',pathname_coreg,'...'])

% CT
    % Load
    CT = niftiread(fullfile([pathname_coreg, scans_coreg(1).name]));
    CT_info = niftiinfo(fullfile([pathname_coreg, scans_coreg(1).name]));
% MRI
    % Load
    MRI = niftiread(fullfile([pathname_coreg, scans_coreg(4).name]));
    MRI_info = niftiinfo(fullfile([pathname_coreg, scans_coreg(4).name]));
% PET
    % Load
    PET = niftiread(fullfile([pathname_coreg, scans_coreg(5).name]));
    PET_info = niftiinfo(fullfile([pathname_coreg, scans_coreg(5).name]));

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IDEAL DOSE MAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load variables
load('BoundingBox.mat')
load('PET_VOI.mat')

% Variables
IDM     = zeros(size(PET));     % Ideal Dose Map (IDM)
I_high  = max(max(max(PET_VOI)));   % TODO maybe better to use maximum in bounding box area, not global max?
I_low   = I_high * 0.25;
D_high  = 2800;                 % unit: cGy
D_low   = 2000;                 % unit: cGy

for i = 1 : size(PET,1)
    for j = 1 : size(PET,2)
        for k = 1 : size(PET,3)
            
            % TODO only do the calculation in bounding box
            Intensity = PET_VOI(i,j,k);
            PresDose = PrescribedDose(Intensity,I_high,I_low,D_high,D_low);
            IDM(i,j,k) = PresDose;
            
        end
    end
end

clearvars i j k Intensity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
orthosliceViewer(IDM)

figure()
orthosliceViewer(PET)

imtool(IDM(:,:,85),[1500 2900])

