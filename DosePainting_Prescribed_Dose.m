%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prescribed Dose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Information
% 
% This script makes a map of the ideal prescribed dose. This dose depends 
% on the signal intensity from the functional imaging. The function below
% indicates the linear relationship that is used to obtain the ideal dose
% map.
%
% D(I) = D(low) + (I - I(low))/(I(high) - I(low))*(D(high) - D(low))
% 
% with
% I             = voxel intensity
% D(I)          = prescribed dose
% D(high)       = 28 Gy
% D(low)        = 20 Gy
% I(high)       = 95th percentile of PET voxel intensity within PTV        
% I(low)        = I(high)*0.25
%
%
% Written by Sam Donche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLEAN SLATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;
clc;
imtool close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS TO ADJUST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Work directory
workdir = pwd;                          % Current MATLAB directory is important!!
% or pathname = 'C:\...';

% Load extra files
boundingbox_file = 'BoundingBox.mat';   % Bounding Box
PET_VOI_file = 'PET_VOI.mat';           % PET VOI
PET_VOI_50_file = 'PET_VOI_50.mat';     % VOI50

% Dose painting values
D_high  = 2800;                         % Highest dose; unit: cGy
D_low   = 2000;                         % Lowest dose; unit: cGy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB TOOLBOXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/download/)
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Extra\spm12')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Preclinical\Scripts\App_BoundingBox')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Preclinical\DosePainting\')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ IMAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Workdirectory
cd(strcat(workdir,'\'));
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
% MRI (coreg)
    % Load
    MRI = niftiread(fullfile([pathname_coreg, scans_coreg(4).name]));
    MRI_info = niftiinfo(fullfile([pathname_coreg, scans_coreg(4).name]));
% PET (coreg)
    % Load
    PET = niftiread(fullfile([pathname_coreg, scans_coreg(5).name]));
    PET_info = niftiinfo(fullfile([pathname_coreg, scans_coreg(5).name]));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IDEAL DOSE MAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load variables
load(boundingbox_file)
load(PET_VOI_file)
load(PET_VOI_50_file)

% Variables
IDM     = zeros(size(PET));                     % Ideal Dose Map (IDM)
I_high  = prctile(nonzeros(PET_VOI_50),95);     % 95 percentile in SUV 50 volume
I_low   = I_high * 0.25;

for i = 1 : size(PET,1)
    for j = 1 : size(PET,2)
        for k = 1 : size(PET,3)
            
            Intensity = PET_VOI_50(i,j,k);
            if Intensity > 0
                PresDose = PrescribedDose(Intensity,I_high,I_low,D_high,D_low);
                IDM(i,j,k) = PresDose;
            end
        end
    end
end

clearvars i j k Intensity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE IDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File for further calculations

save('IDM.mat','IDM')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WRITE NIFTI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File to view the ideal dose map in external software

IDM16 = int16(IDM);
niftiwrite(IDM16,'IDM.nii',CT_info)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VIEW IMAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    figure()
    orthosliceViewer(IDM)

    figure()
    orthosliceViewer(PET)
end

