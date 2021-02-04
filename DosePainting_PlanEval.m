%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plan Evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Information
%
% This script generates a dose-volume histogram (DVH), Q-volume histogram
% (QVH) and the Q-factor. These parameters can be used to evaluate the
% radiation plan. 
%
% A quality index for each pixel Q(p) is obtained by dividing the obtained
% dose over the intended dose. These values can be visualised in a Q-volume
% histogram, by displaying the partial volume for which the QF is greater
% than or equal to each avscis value. Ideally such a curve would drop
% steeply at Q = 1. The quality-factor (QF) is defined as the mean absolute
% deviation of Q to 1 within the given volume:
%
% QF = (1/n)*sum(p)abs(Q(p)-1)
%
% This method is adopted from [18F]FDG PET voxel intensity-based IMRT for
% head and neck cancer at Ghent University Hospital
% (https://doi.org/10.1016/j.radonc.2006.03.003)
%
% Written by Sam Donche

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
workdir = pwd;                  % Current MATLAB directory is important!!
% or pathname = 'C:\...';

% Directory with NIfTI file SARRP dose, called "SARRPDose.nii"
dosedir = 'Method8+1';

% Load extra files
IDM_file = 'IDM.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB TOOLBOXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/download/)
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
% MRI
    % Load
    MRI = niftiread(fullfile([pathname_coreg, scans_coreg(4).name]));
    MRI_info = niftiinfo(fullfile([pathname_coreg, scans_coreg(4).name]));
% PET
    % Load
    PET = niftiread(fullfile([pathname_coreg, scans_coreg(5).name]));
    PET_info = niftiinfo(fullfile([pathname_coreg, scans_coreg(5).name]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ DOSE MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load SARRP dose map
file_name_map = strcat(dosedir,'\SARRPDose.nii');
X = ['Reading ', file_name_map,'...'];
disp(X); clearvars X;
DoseMap = niftiread(fullfile([pathname,file_name_map]));
DoseMap_info = niftiinfo(fullfile([pathname,file_name_map]));

% Load ideal dose map
load(IDM_file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REORIENT DOSEMAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reorient dose map calculated by SARRP

DoseMap_flipped = flip(DoseMap,2);

if 0
    figure()
    orthosliceViewer(DoseMap_flipped)
    figure()
    orthosliceViewer(DoseMap)
    figure()
    orthosliceViewer(IDM)
    figure()
    orthosliceViewer(CT)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE Q-factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qp = zeros(nnz(IDM),1);
counter = 1;

for i = 1 : size(PET,1)
    for j = 1 : size(PET,2)
        for k = 1 : size(PET,3)
        
            IdealDose = IDM(i,j,k);
            if IdealDose > 1
                
                ObtainedDose = DoseMap_flipped(i,j,k);
                Qp(counter) = ObtainedDose/IdealDose;
                counter = counter + 1;
                
            end
        end
    end
end
    
Q_factor = sum(abs(Qp -1))/length(Qp);
display(Q_factor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Q-VOLUME HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q_hist = 0:0.001:1.5;
Vol_hist = zeros(1,length(Q_hist));

for i = 1 : length(Q_hist)
    Q_val = Q_hist(i);
    Vol_hist(i) = sum(Qp >= Q_val)/length(Qp)*100;    
end

hist = scatter(Q_hist,Vol_hist);
xlabel('Q-value')
ylabel('Volume (%)')

saveas(hist,strcat(dosedir,'\Q-Volume_Histogram.fig'))
saveas(hist,strcat(dosedir,'\Q-Volume_Histogram.jpg'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DIFFERENCE MAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra info
%       DoseMap too high --> negative value
%       DoseMap too low  --> positive value

Diff = IDM - DoseMap_flipped;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    figure()
    orthosliceViewer(DoseMap_flipped)
    colormap(jet)

    figure()
    orthosliceViewer(IDM)
    colormap(jet)

    figure()
    orthosliceViewer(Diff)
    colormap(jet)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE NIFTI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Files to view the ideal dose map in external software

% Dose Map
niftiwrite(DoseMap_flipped,strcat('DoseMap','_',dosedir,'.nii'),DoseMap_info)

% Difference Map
Diff16 = int16(Diff);
niftiwrite(Diff16,strcat('DifferenceMap','_',dosedir,'.nii'),CT_info)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOSE-VOLUME HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Doses = zeros(nnz(IDM),1);
counter = 1;

for i = 1 : size(PET,1)
    for j = 1 : size(PET,2)
        for k = 1 : size(PET,3)
        
            IdealDose = IDM(i,j,k);
            if IdealDose > 1
                
                Doses(counter) = DoseMap_flipped(i,j,k);
                Ideal_Doses(counter) = IdealDose;
                counter = counter + 1;
                
            end
        end
    end
end

D_hist = 0:0.01:3500;
DVol_hist = zeros(1,length(D_hist));

for i = 1 : length(D_hist)
    
    Dose = D_hist(i);
    DVol_hist(i) = sum(Doses >= Dose)/length(Doses)*100;
    
end

D_hist_IDM = 0:0.01:3500;
DVol_hist_IDM = zeros(1,length(D_hist));

for i = 1 : length(D_hist_IDM)
    
    Dose = D_hist_IDM(i);
    DVol_hist_IDM(i) = sum(Ideal_Doses >= Dose)/length(Ideal_Doses)*100;
    
end

hist2 = scatter(D_hist,DVol_hist);
hold on
scatter(D_hist_IDM,DVol_hist_IDM);
xlabel('Dose')
ylabel('Volume (%)')

saveas(hist2,strcat(dosedir,'\Dose-Volume_Histogram.fig'))
saveas(hist2,strcat(dosedir,'\Dose-Volume_Histogram.jpg'))

