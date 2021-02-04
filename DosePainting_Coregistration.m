%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Dose Painting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This script coregisters three imaging modalities (CT, MRI and PET).
% It consists of multiple steps:
%       - Conversion of DICOM to NIfTI
%       - Filter image (PET)
%       - Adjust orientation
%       - Image cropping (PET)
%       - Move image center (PET)
%       - Coregistration
%
% Script usage:
%   - Move MATLAB to a workdirectory
%   - Move images to seperate folders in the workdirectory
%   - Name the folders "CT", "MRI" and "PET"
%
%   Workdirectory
%       -> CT
%       -> MRI
%       -> PET
%
% Animal orientation during imaging
%   - CBCT from SARRP (XStrahl) = REFERENCE IMAGE
%       -> head towards entrance, move motorized bed 10 mm downwards
%   - MRI from Biospin Pharmascan (Bruker)
%       -> head first prone
%   - PET from Beta-cube (Molecubes)
%       -> default
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
%% MATLAB TOOLBOXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/download/)
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Extra\spm12')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Preclinical\Scripts\App_BoundingBox')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Preclinical\Scripts')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VARIABLES TO ADJUST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Workdirectory
pathname = pwd;
% or pathname = 'C:\...';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONVERSION DICOM TO NIFTI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For conversion from DICOM to NIfTI use the function dcm2niix from mricron
% MATLAB can be used to perform tasks on command line level
% Use "!" before the command 
% help: !dcm2niix -h &

tic                                 % START TIMER

% Make directory
!mkdir 1NIfTI
disp('Making NIfTI directory')

% Convert CT
    disp('Converting CT')
    !dcm2niix -o "./1NIfTI/" -c "CT" -d 0 -s y -f 1CT_%i ./CT
    disp(' ')

% Convert MRI
    disp('Converting MRI')
    !dcm2niix -o "./1NIfTI/" -c "MR" -d 0 -s y -f 2MRI_%n ./MRI
    disp(' ')

% Convert PET
    disp('Converting PET')
    !dcm2niix -o "./1NIfTI/" -c "PT" -d 0 -s y -f 3PET_%i ./PET
    disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ IMAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc                         % clear command window

%Data directory
pathname_nii = [pathname,'\','1NIfTI\'];
nii = dir(pathname_nii);
nii = nii(3:end);           % Remove . and ..

% Load with matlab functions
disp(['Reading files from ',pathname_nii,'...'])

% CT
    % Load
    CT_json = jsondecode(fileread(fullfile([pathname_nii, nii(1).name])));
    CT = niftiread(fullfile([pathname_nii, nii(2).name]));
    CT_info = niftiinfo(fullfile([pathname_nii, nii(2).name]));

    % Check
    if strcmp('CT',CT_json.Modality) == 0 || strcmp('CT',CT_info.raw.aux_file) == 0
        warning('CT and CT_info variables do not contain a CT image')
    end

% MRI
    % Load
    MRI_json = jsondecode(fileread(fullfile([pathname_nii, nii(3).name])));
    MRI = niftiread(fullfile([pathname_nii, nii(4).name]));
    MRI_info = niftiinfo(fullfile([pathname_nii, nii(4).name]));

    % Check
    if strcmp('MR',MRI_json.Modality) == 0 || strcmp('MR',MRI_info.raw.aux_file) == 0
        warning('MRI and MRI_info variables do not contain a MR image')
    end

% PET
    % Load
    PET_json = jsondecode(fileread(fullfile([pathname_nii, nii(5).name])));
    PET = niftiread(fullfile([pathname_nii, nii(6).name]));
    PET_info = niftiinfo(fullfile([pathname_nii, nii(6).name]));

    % Check
    if strcmp('PT',PET_json.Modality) == 0 || strcmp('PT',PET_info.raw.aux_file) == 0
        warning('PET and PET_info variables do not contain a PET image')
    end

% Qform is deleted from header TODO: write function for this
% CT
    CT_info.raw.qform_code = 0;
    CT_info.raw.quatern_b = 0;
    CT_info.raw.quatern_c = 0;
    CT_info.raw.quatern_d = 0;
    CT_info.raw.qoffset_x = 0;
    CT_info.raw.qoffset_y = 0;
    CT_info.raw.qoffset_z = 0;

% MRI
    MRI_info.raw.qform_code = 0;
    MRI_info.raw.quatern_b = 0;
    MRI_info.raw.quatern_c = 0;
    MRI_info.raw.quatern_d = 0;
    MRI_info.raw.qoffset_x = 0;
    MRI_info.raw.qoffset_y = 0;
    MRI_info.raw.qoffset_z = 0;
    
% PET
    PET_info.raw.qform_code = 0;
    PET_info.raw.quatern_b = 0;
    PET_info.raw.quatern_c = 0;
    PET_info.raw.quatern_d = 0;
    PET_info.raw.qoffset_x = 0;
    PET_info.raw.qoffset_y = 0;
    PET_info.raw.qoffset_z = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FILTER PET DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PET data reconstruction:
% OSEM 30 iterations
%
% Amide filter: 
% Kernel size = 31 (pixels?)
% FWHM(mm) = 1 mm

FWHM = 1;
sigma = FWHM / (sqrt(8 * log(2))* PET_info.PixelDimensions(1));
PET = imgaussfilt3(PET, sigma, 'FilterSize', 31);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADJUST ORIENTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make directory
!mkdir 2Reorient

% CT
    % Flip around Z-axis
    CT_info.Transform.T(3,1:3) = -(CT_info.Transform.T(3,1:3));
    CT_info.Transform.T(4,2) = -(CT_info.Transform.T(4,2));

    % Save new nifti file
    niftiwrite(CT,'2Reorient\1CT_reorient.nii',CT_info);

% MRI
    % Flip around X-axis
    MRI_info.Transform.T(1,1:3) = -(MRI_info.Transform.T(1,1:3));
    MRI_info.Transform.T(4,1) = -(MRI_info.Transform.T(4,1));
 
    % Save new nifti file
    niftiwrite(MRI,'2Reorient\2MRI_reorient.nii',MRI_info)

% PET
    % Flip around X- and Y-axis
    PET_info.Transform.T(2,1:3) = -(PET_info.Transform.T(2,1:3));
    PET_info.Transform.T(4,2) = -(PET_info.Transform.T(4,2));
     
    % Save new nifti file
    niftiwrite(PET,'2Reorient\3PET_reorient.nii',PET_info)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMAGE CROPPING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear previous data
clearvars -except pathname pathname_coreg

% Path
pathname_reo = [pathname,'\2Reorient\'];
scans_reo = dir(pathname_reo);
scans_reo = scans_reo(3:end); % Remove . and ..

% Make directory and copy non-cropped images
!mkdir 3Crop
!copy /Y 2Reorient\1CT_reorient.nii "3Crop\1CT_noncropped.nii"
!copy /Y 2Reorient\2MRI_reorient.nii "3Crop\2MRI_noncropped.nii"

% Parameters
% PET
% # pixels removed from each side
crop_pet_x_1 = 40;      % Left side
crop_pet_x_2 = 40;      % Right side
crop_pet_y_1 = 60;      % Dorsaal side
crop_pet_y_2 = 40;      % Ventral side
crop_pet_z_1 = 170;     % Inferior
crop_pet_z_2 = 30;      % Superior

% Load PET
PET = niftiread(fullfile([pathname_reo, scans_reo(3).name]));
PET_info = niftiinfo(fullfile([pathname_reo, scans_reo(3).name]));

% Crop PET
[pet_x, pet_y, pet_z] = size(PET);

    % Change data matrix
    % X
    PET(pet_x - crop_pet_x_1 + 1 : pet_x,:,:) = [];
    PET(1:crop_pet_x_2,:,:) = [];
    
    % Y
    PET(:,pet_y - crop_pet_y_1 + 1 : pet_y,:) = [];
    PET(:,1:crop_pet_y_2,:) = [];
   
    % Z
    PET(:,:,pet_z - crop_pet_z_1 + 1 : pet_z) = [];
    PET(:,:,1:crop_pet_z_2) = [];

    % Adjust header
    PET_info.ImageSize(1) = pet_x-crop_pet_x_1-crop_pet_x_2;
    PET_info.ImageSize(2) = pet_y-crop_pet_y_1-crop_pet_y_2;
    PET_info.ImageSize(3) = pet_z-crop_pet_z_1-crop_pet_z_2;

    PET_info.raw.dim(2) = pet_x-crop_pet_x_1-crop_pet_x_2;
    PET_info.raw.dim(3) = pet_y-crop_pet_y_1-crop_pet_y_2;
    PET_info.raw.dim(4) = pet_z-crop_pet_z_1-crop_pet_z_2;

% Save
niftiwrite(PET,'3Crop\3PET_cropped.nii',PET_info)    
    
% Check crop
    figure()
    orthosliceViewer(PET)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MOVE IMAGE CENTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make directory and copy non-cropped images
!mkdir 4Trans
!copy /Y 3Crop\1CT_noncropped.nii "4Trans\1CT_nontrans.nii"
!copy /Y 3Crop\2MRI_noncropped.nii "4Trans\2MRI_nontrans.nii"

% Path
pathname_crop = [pathname,'\3Crop\'];
scans_crop = dir(pathname_crop);
scans_crop = scans_crop(3:end); % Remove . and ..

% CT
    % Load
    CT = niftiread(fullfile([pathname_crop, scans_crop(1).name]));
    CT_info = niftiinfo(fullfile([pathname_crop, scans_crop(1).name]));

    % Extract info
    pix_ct_x = sum(CT_info.Transform.T(1,1:3));
    pix_ct_z = sum(CT_info.Transform.T(2,1:3));
    pix_ct_y = sum(CT_info.Transform.T(3,1:3));
    dim_ct_x = CT_info.ImageSize(1);
    dim_ct_z = CT_info.ImageSize(2);
    dim_ct_y = CT_info.ImageSize(3);
    delta_ct_x = CT_info.Transform.T(4,1);
    delta_ct_z = CT_info.Transform.T(4,2);
    delta_ct_y = CT_info.Transform.T(4,3);
    
% PET
    % Load
    PET = niftiread(fullfile([pathname_crop, scans_crop(3).name]));
    PET_info = niftiinfo(fullfile([pathname_crop, scans_crop(3).name]));
    
    % Extract info
    pix_pet_x = sum(PET_info.Transform.T(1,1:3));
    pix_pet_y = sum(PET_info.Transform.T(2,1:3));
    pix_pet_z = sum(PET_info.Transform.T(3,1:3));
    dim_pet_x = PET_info.ImageSize(1);
    dim_pet_y = PET_info.ImageSize(2);
    dim_pet_z = PET_info.ImageSize(3);    
    
    % Adjust PET translation values; TODO: aanpassen wegens niet algemeen
    delta_pet_x = abs(dim_pet_x*pix_pet_x/2);
    delta_pet_y = abs(dim_ct_y*pix_ct_y*2/5)-abs(dim_ct_y*pix_ct_y/2);
    delta_pet_z = abs(dim_pet_z*pix_pet_z/2);
    
    PET_info.Transform.T(4,1) = delta_pet_x;
    PET_info.Transform.T(4,2) = delta_pet_y;
    PET_info.Transform.T(4,3) = delta_pet_z;
    
    PET_info.raw.srow_x(4) = delta_pet_x;
    PET_info.raw.srow_y(4) = delta_pet_y;
    PET_info.raw.srow_z(4) = delta_pet_z;
    
    niftiwrite(PET,'4Trans\3PET_trans.nii',PET_info)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMAGE REGISTRATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TODO: original files (e.g. 2MRI of 3PET) worden ook aangepast! waarom?

% Clear workspace
clearvars -except pathname pathname_reo

% Copy directory
    !if exist 5Coregister rmdir /S /Q 5Coregister & echo 5Coregister has been deleted
    !mkdir "5Coregister"
    !copy /Y 4Trans "5Coregister"
    
% Path
%Data directory
pathname_coreg = [pathname,'\','5Coregister\'];
scans_coreg = dir(pathname_coreg);
scans_coreg = scans_coreg(3:end); % Remove . and ..

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%Coregistratie CT-MRI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global ref_image warp_image

    ref_image = [pathname_coreg,scans_coreg(1).name,',1']; %CT
    warp_image = [pathname_coreg,scans_coreg(2).name, ',1']; %MRI

    % List of open inputs
    nrun = 5; % enter the number of runs here 
    jobfile = {'C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Preclinical\DosePainting\coregister_preclinical_job.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(0, nrun);
    for crun = 1:nrun
    end
    spm('defaults', 'PET');
    spm_jobman('run', jobs, inputs{:});
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Coregistratie CT - PET%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    warp_image = [pathname_coreg,scans_coreg(3).name, ',1']; %PET

    % List of open inputs
    nrun = 10; % enter the number of runs here
    jobfile = {'C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Preclinical\DosePainting\coregister_preclinical_job_PET.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(0, nrun);
    for crun = 1:nrun
    end
    spm('defaults', 'PET');
    spm_jobman('run', jobs, inputs{:});

    toc                         % STOP TIMER
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Examine coregistration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Often better to check with external imaging software (e.g. Amide)
% And adjust coregistration if necessary

% Clear workspace
clearvars -except pathname pathname_coreg  

% Path
scans_coreg = dir(pathname_coreg);
scans_coreg = scans_coreg(3:end); % Remove . and ..

% Load
CT = niftiread(fullfile([pathname_coreg, scans_coreg(1).name]));
CT_info = niftiinfo(fullfile([pathname_coreg, scans_coreg(1).name]));  

MRI_coreg = niftiread(fullfile([pathname_coreg, scans_coreg(4).name]));
MRI_coreg_info = niftiinfo(fullfile([pathname_coreg, scans_coreg(4).name]));

PET_coreg = niftiread(fullfile([pathname_coreg, scans_coreg(5).name]));
PET_coreg_info = niftiinfo(fullfile([pathname_coreg, scans_coreg(5).name]));

% Display

CmpCT_MRI = CompImages(CT,MRI_coreg);

waitfor(CmpCT_MRI)

CmpCT_PET = CompImages(CT,PET_coreg);
