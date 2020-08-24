%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Dose Painting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Information
%
%        Maximum field size: 40 x 80 mm 
%        Minimum field size: 1 x 1 mm
%        SARRP resolution: 0.01 mm
%
% Conform PhD Sam Donche

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLEAN SLATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;
clc;
imtool close all;

% TODO bed 0° en andere mogelijkheden
% TODO literatuur variabele collimator
% TODO methode voor dosis nagaan
% TODO vergelijking MRI
% TODO make json file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB TOOLBOXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/download/)
%   NIFTI and ANALYZE tools (https://nl.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
%   Medical Image Reader and Viewer (https://www.mathworks.com/matlabcentral/fileexchange/53745-medical-image-reader-and-viewer)
%   nifti_utils (https://github.com/justinblaber/nifti_utils)

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
%% DEFINE VOI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cd 'C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Preclinical\Scripts'

DP_output = DosePainting;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE CONTOURS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract VOI from PET
    PET_VOI = double(PET) .* boundingbox_var;
    save('PET_VOI.mat','PET_VOI')
    save('BoundingBox.mat','boundingbox_var')
    %figure(); sliceViewer(PET)
    %figure(); sliceViewer(boundingbox_var)
    %figure(); sliceViewer(PET_VOI)
    figure(); h = imagesc(squeeze(PET_VOI(:,:,85))); xlabel('X-as'); ylabel('Y-as'); colormap gray; 
    hold on
    plot(240,207,'or')

% PET max in VOI
    PET_VOI_MAX = max(max(max(PET_VOI)));

% Calculate VOI 50
    % PET 50 threshold
        THRES_PET_50 = PET_VOI_MAX*0.50;

    % PET 50 thresholds    
        PET_VOI_SEG_50 = imbinarize(PET_VOI,THRES_PET_50);
        %figure(); sliceViewer(PET_VOI_SEG_50)

    % FILL holes VOI 50
        PET_VOI_SEG_FILL_50 = imfill(PET_VOI_SEG_50, 'holes');

    % Extract largest segment
        PET_VOI_SEG_LARG_50 = ExtractSegments(PET_VOI_SEG_FILL_50,1);
    
% Use VOI 50 as a 'bounding box' for other VOI definitions
    PET_VOI_50 = double(PET) .* PET_VOI_SEG_LARG_50;
    save('PET_VOI_50.mat','PET_VOI_50')

% Calculate other VOIs
% PET thresholds
    THRES_PET_60 = PET_VOI_MAX*0.60;
    THRES_PET_70 = PET_VOI_MAX*0.70;
    THRES_PET_80 = PET_VOI_MAX*0.80;
    THRES_PET_90 = PET_VOI_MAX*0.90;
    THRES_PET_95 = PET_VOI_MAX*0.95;

% PET segmentations
    PET_VOI_SEG_60 = imbinarize(PET_VOI_50,THRES_PET_60);
    PET_VOI_SEG_70 = imbinarize(PET_VOI_50,THRES_PET_70);
    PET_VOI_SEG_80 = imbinarize(PET_VOI_50,THRES_PET_80);
    PET_VOI_SEG_90 = imbinarize(PET_VOI_50,THRES_PET_90);
    PET_VOI_SEG_95 = imbinarize(PET_VOI_50,THRES_PET_95);

% Fill holes
    PET_VOI_SEG_FILL_50 = imfill(PET_VOI_SEG_50, 'holes');
    PET_VOI_SEG_FILL_60 = imfill(PET_VOI_SEG_60, 'holes');
    PET_VOI_SEG_FILL_70 = imfill(PET_VOI_SEG_70, 'holes');
    PET_VOI_SEG_FILL_80 = imfill(PET_VOI_SEG_80, 'holes');
    PET_VOI_SEG_FILL_90 = imfill(PET_VOI_SEG_90, 'holes');
    PET_VOI_SEG_FILL_95 = imfill(PET_VOI_SEG_95, 'holes');

% Display result
    if 0
        cmpimg1 = CompImages(PET_VOI,PET_VOI_SEG_LARG_50);
        waitfor(cmpimg1);
        cmpimg2 = CompImages(PET_VOI,PET_VOI_SEG_FILL_60);
        waitfor(cmpimg2);
        cmpimg3 = CompImages(PET_VOI,PET_VOI_SEG_FILL_70);
        waitfor(cmpimg3);
        cmpimg4 = CompImages(PET_VOI,PET_VOI_SEG_FILL_80);
        waitfor(cmpimg4);
        cmpimg5 = CompImages(PET_VOI,PET_VOI_SEG_FILL_90);
        waitfor(cmpimg5);
        cmpimg6 = CompImages(PET_VOI,PET_VOI_SEG_FILL_95);
        waitfor(cmpimg6);
        clear cmpimg1 cmpimg2 cmpimg3 cmpimg4 cmpimg5 cmpimg6
    end

% Information on VOIs
    VOI = [];
    [labeledImage_50, numberOfBlobs_50] = bwlabeln(PET_VOI_SEG_LARG_50);
    blobMeasurements_PET_50 = regionprops3(labeledImage_50, PET_VOI, 'Volume', 'Centroid', 'WeightedCentroid','MaxIntensity','MeanIntensity','MinIntensity','Orientation','BoundingBox','EquivDiameter');
    for i = 1:numberOfBlobs_50
        VOI = [VOI;50];
    end
    blobMeasurements_PET_50 = addvars(blobMeasurements_PET_50,VOI,'Before','Volume');
    
    VOI = [];
    [labeledImage_60, numberOfBlobs_60] = bwlabeln(PET_VOI_SEG_FILL_60);
    blobMeasurements_PET_60 = regionprops3(labeledImage_60, PET_VOI, 'Volume', 'Centroid', 'WeightedCentroid','MaxIntensity','MeanIntensity','MinIntensity','Orientation','BoundingBox','EquivDiameter');
    for i = 1:numberOfBlobs_60
        VOI = [VOI;60];
    end
    blobMeasurements_PET_60 = addvars(blobMeasurements_PET_60,VOI,'Before','Volume');
    
    VOI = [];
    [labeledImage_70, numberOfBlobs_70] = bwlabeln(PET_VOI_SEG_FILL_70);
    blobMeasurements_PET_70 = regionprops3(labeledImage_70, PET_VOI, 'Volume', 'Centroid', 'WeightedCentroid','MaxIntensity','MeanIntensity','MinIntensity','Orientation','BoundingBox','EquivDiameter');
    for i = 1:numberOfBlobs_70
        VOI = [VOI;70];
    end
    blobMeasurements_PET_70 = addvars(blobMeasurements_PET_70,VOI,'Before','Volume');
    
    
    VOI = [];
    [labeledImage_80, numberOfBlobs_80] = bwlabeln(PET_VOI_SEG_FILL_80);
    blobMeasurements_PET_80 = regionprops3(labeledImage_80, PET_VOI, 'Volume', 'Centroid', 'WeightedCentroid','MaxIntensity','MeanIntensity','MinIntensity','Orientation','BoundingBox','EquivDiameter');
    for i = 1:numberOfBlobs_80
        VOI = [VOI;80];
    end
    blobMeasurements_PET_80 = addvars(blobMeasurements_PET_80,VOI,'Before','Volume');
    
    VOI = [];
    [labeledImage_90, numberOfBlobs_90] = bwlabeln(PET_VOI_SEG_FILL_90);
    blobMeasurements_PET_90 = regionprops3(labeledImage_90, PET_VOI, 'Volume', 'Centroid', 'WeightedCentroid','MaxIntensity','MeanIntensity','MinIntensity','Orientation','BoundingBox','EquivDiameter');
    for i = 1:numberOfBlobs_90
        VOI = [VOI;90];
    end
    blobMeasurements_PET_90 = addvars(blobMeasurements_PET_90,VOI,'Before','Volume');
    
    VOI = [];
    [labeledImage_95, numberOfBlobs_95] = bwlabeln(PET_VOI_SEG_FILL_95);
    blobMeasurements_PET_95 = regionprops3(labeledImage_95, PET_VOI, 'Volume', 'Centroid', 'WeightedCentroid','MaxIntensity','MeanIntensity','MinIntensity','Orientation','BoundingBox','EquivDiameter');
    for i = 1:numberOfBlobs_95
        VOI = [VOI;95];
    end
    blobMeasurements_PET_95 = addvars(blobMeasurements_PET_95,VOI,'Before','Volume');
    
    vertcat(blobMeasurements_PET_50,blobMeasurements_PET_60,blobMeasurements_PET_70,blobMeasurements_PET_80,...
        blobMeasurements_PET_90,blobMeasurements_PET_95)

    clearvars labeledImage_50 labeledImage_60 labeledImage_70 labeledImage_80 labeledImage_90 labeledImage_95 ...
        THRES_PET_50 THRES_PET_60 THRES_PET_70 THRES_PET_80 THRES_PET_90 THRES_PET_95 ...
        PET_VOI_SEG_50 PET_VOI_SEG_60 PET_VOI_SEG_70 PET_VOI_SEG_80 PET_VOI_SEG_90 PET_VOI_SEG_95 ...
        PET_VOI_MAX VOI
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE ISOCENTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate center of multiple blobs if necessary
    
    % Extract centroids
    VOI50_centroid = blobMeasurements_PET_50.Centroid;
    VOI60_centroid = blobMeasurements_PET_60.Centroid;
    VOI70_centroid = blobMeasurements_PET_70.Centroid;
    VOI80_centroid = blobMeasurements_PET_80.Centroid;
    VOI90_centroid = blobMeasurements_PET_90.Centroid;
    VOI95_centroid = blobMeasurements_PET_95.Centroid;

    % Change centroid
    if numberOfBlobs_60 > 1
        VOI60_centroid = mean(blobMeasurements_PET_60.Centroid,1);
    end
    if numberOfBlobs_70 > 1
        VOI70_centroid = mean(blobMeasurements_PET_70.Centroid,1);
    end
    if numberOfBlobs_80 > 1
        VOI80_centroid = mean(blobMeasurements_PET_80.Centroid,1);
    end
    if numberOfBlobs_90 > 1
        VOI90_centroid = mean(blobMeasurements_PET_90.Centroid,1);
    end
    if numberOfBlobs_95 > 1
        VOI95_centroid = mean(blobMeasurements_PET_95.Centroid,1);
    end

    clearvars numberOfBlobs_50 numberOfBlobs_60 numberOfBlobs_70 numberOfBlobs_80 numberOfBlobs_90 numberOfBlobs_95
        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROJECT VOIs ON PERPENDICULAR PLANE TO INCIDENT BEAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Planes are defined by a normal vector and a point.
% The normal vector defined by the gantry angle
% The point is defined by the centroid of the VOI
% Equation plane:
% a*(x-x1) + b*(y-y1) + c*(z-z1) = 0
% with
% (a,b,c) the normal vector
% (x1,y1,z1) the point

% Normal vector to construct plane
norm_0 = [0,0,1];               % Gantry 0°
norm_40 = [1,0,1/tand(40)];     % Gantry 40°
norm_80 = [1,0,1/tand(80)];     % Gantry 80°
norm_120 = [1,0,-tand(30)];     % Gantry 120°

% Original points
    % VOI50
    % TODO make size variables to use instead of using different once all
    % the time
    ori_index_50 = find(PET_VOI_SEG_LARG_50);
    [original_y_50, original_x_50, original_z_50] = ind2sub((size(PET_VOI_SEG_LARG_50)), ori_index_50); % ind2sub output = [ROW, COLUMN,]
    
    % VOI60
    ori_index_60 = find(PET_VOI_SEG_FILL_60);
    [original_y_60, original_x_60, original_z_60] = ind2sub((size(PET_VOI_SEG_FILL_60)), ori_index_60); % ind2sub output = [ROW, COLUMN,]
    
    % VOI70
    ori_index_70 = find(PET_VOI_SEG_FILL_70);
    [original_y_70, original_x_70, original_z_70] = ind2sub((size(PET_VOI_SEG_FILL_70)), ori_index_70); % ind2sub output = [ROW, COLUMN,]
    
    % VOI80
    ori_index_80 = find(PET_VOI_SEG_FILL_80);
    [original_y_80, original_x_80, original_z_80] = ind2sub((size(PET_VOI_SEG_FILL_80)), ori_index_80); % ind2sub output = [ROW, COLUMN,]
    
    % VOI90
    ori_index_90 = find(PET_VOI_SEG_FILL_90);
    [original_y_90, original_x_90, original_z_90] = ind2sub((size(PET_VOI_SEG_FILL_90)), ori_index_90); % ind2sub output = [ROW, COLUMN,]
    
    % VOI95
    ori_index_95 = find(PET_VOI_SEG_FILL_95);
    [original_y_95, original_x_95, original_z_95] = ind2sub((size(PET_VOI_SEG_FILL_95)), ori_index_95); % ind2sub output = [ROW, COLUMN,]
    
    clearvars ori_index_50 ori_index_60 ori_index_70 ori_index_80 ori_index_90 ori_index_95

% Projection
    % Declaration
    % VOI50
        projection_50_n0_x = zeros(length(original_x_50),1);
        projection_50_n0_y = zeros(length(original_y_50),1);
        projection_50_n0_z = zeros(length(original_z_50),1);
        
        projection_50_n40_x = zeros(length(original_x_50),1);
        projection_50_n40_y = zeros(length(original_y_50),1);
        projection_50_n40_z = zeros(length(original_z_50),1);
        
        projection_50_n80_x = zeros(length(original_x_50),1);
        projection_50_n80_y = zeros(length(original_y_50),1);
        projection_50_n80_z = zeros(length(original_z_50),1);
        
        projection_50_n120_x = zeros(length(original_x_50),1);
        projection_50_n120_y = zeros(length(original_y_50),1);
        projection_50_n120_z = zeros(length(original_z_50),1);
        
    % VOI60
        projection_60_n0_x = zeros(length(original_x_60),1);
        projection_60_n0_y = zeros(length(original_y_60),1);
        projection_60_n0_z = zeros(length(original_z_60),1);
        
        projection_60_n40_x = zeros(length(original_x_60),1);
        projection_60_n40_y = zeros(length(original_y_60),1);
        projection_60_n40_z = zeros(length(original_z_60),1);
        
        projection_60_n80_x = zeros(length(original_x_60),1);
        projection_60_n80_y = zeros(length(original_y_60),1);
        projection_60_n80_z = zeros(length(original_z_60),1);
        
        projection_60_n120_x = zeros(length(original_x_60),1);
        projection_60_n120_y = zeros(length(original_y_60),1);
        projection_60_n120_z = zeros(length(original_z_60),1);
    
    % VOI70
        projection_70_n0_x = zeros(length(original_x_70),1);
        projection_70_n0_y = zeros(length(original_y_70),1);
        projection_70_n0_z = zeros(length(original_z_70),1);
        
        projection_70_n40_x = zeros(length(original_x_70),1);
        projection_70_n40_y = zeros(length(original_y_70),1);
        projection_70_n40_z = zeros(length(original_z_70),1);
        
        projection_70_n80_x = zeros(length(original_x_70),1);
        projection_70_n80_y = zeros(length(original_y_70),1);
        projection_70_n80_z = zeros(length(original_z_70),1);
        
        projection_70_n120_x = zeros(length(original_x_70),1);
        projection_70_n120_y = zeros(length(original_y_70),1);
        projection_70_n120_z = zeros(length(original_z_70),1);        
    
    % VOI80
        projection_80_n0_x = zeros(length(original_x_80),1);
        projection_80_n0_y = zeros(length(original_y_80),1);
        projection_80_n0_z = zeros(length(original_z_80),1);
        
        projection_80_n40_x = zeros(length(original_x_80),1);
        projection_80_n40_y = zeros(length(original_y_80),1);
        projection_80_n40_z = zeros(length(original_z_80),1);
        
        projection_80_n80_x = zeros(length(original_x_80),1);
        projection_80_n80_y = zeros(length(original_y_80),1);
        projection_80_n80_z = zeros(length(original_z_80),1);
        
        projection_80_n120_x = zeros(length(original_x_80),1);
        projection_80_n120_y = zeros(length(original_y_80),1);
        projection_80_n120_z = zeros(length(original_z_80),1);
        
    % VOI90
        projection_90_n0_x = zeros(length(original_x_90),1);
        projection_90_n0_y = zeros(length(original_y_90),1);
        projection_90_n0_z = zeros(length(original_z_90),1);
        
        projection_90_n40_x = zeros(length(original_x_90),1);
        projection_90_n40_y = zeros(length(original_y_90),1);
        projection_90_n40_z = zeros(length(original_z_90),1);
        
        projection_90_n80_x = zeros(length(original_x_90),1);
        projection_90_n80_y = zeros(length(original_y_90),1);
        projection_90_n80_z = zeros(length(original_z_90),1);
        
        projection_90_n120_x = zeros(length(original_x_90),1);
        projection_90_n120_y = zeros(length(original_y_90),1);
        projection_90_n120_z = zeros(length(original_z_90),1);
    
    % VOI95
        projection_95_n0_x = zeros(length(original_x_95),1);
        projection_95_n0_y = zeros(length(original_y_95),1);
        projection_95_n0_z = zeros(length(original_z_95),1);
        
        projection_95_n40_x = zeros(length(original_x_95),1);
        projection_95_n40_y = zeros(length(original_y_95),1);
        projection_95_n40_z = zeros(length(original_z_95),1);
        
        projection_95_n80_x = zeros(length(original_x_95),1);
        projection_95_n80_y = zeros(length(original_y_95),1);
        projection_95_n80_z = zeros(length(original_z_95),1);
        
        projection_95_n120_x = zeros(length(original_x_95),1);
        projection_95_n120_y = zeros(length(original_y_95),1);
        projection_95_n120_z = zeros(length(original_z_95),1);
        
    % Projection 
    % VOI50
    for i = 1:length(original_x_50)
        
        [p_n0_x, p_n0_y, p_n0_z] = Orthproject(norm_0(1),norm_0(2),norm_0(3),...
            VOI50_centroid(1),VOI50_centroid(2),VOI50_centroid(3),...
            original_x_50(i),original_y_50(i),original_z_50(i));
        
        [p_n40_x, p_n40_y, p_n40_z] = Orthproject(norm_40(1),norm_40(2),norm_40(3),...
            VOI50_centroid(1),VOI50_centroid(2),VOI50_centroid(3),...
            original_x_50(i),original_y_50(i),original_z_50(i));
        
        [p_n80_x, p_n80_y, p_n80_z] = Orthproject(norm_80(1),norm_80(2),norm_80(3),...
            VOI50_centroid(1),VOI50_centroid(2),VOI50_centroid(3),...
            original_x_50(i),original_y_50(i),original_z_50(i));
        
        [p_n120_x, p_n120_y, p_n120_z] = Orthproject(norm_120(1),norm_120(2),norm_120(3),...
            VOI50_centroid(1),VOI50_centroid(2),VOI50_centroid(3),...
            original_x_50(i),original_y_50(i),original_z_50(i));
        
        projection_50_n0_x(i) = p_n0_x;
        projection_50_n0_y(i) = p_n0_y;
        projection_50_n0_z(i) = p_n0_z;
        
        projection_50_n40_x(i) = p_n40_x;
        projection_50_n40_y(i) = p_n40_y;
        projection_50_n40_z(i) = p_n40_z;
        
        projection_50_n80_x(i) = p_n80_x;
        projection_50_n80_y(i) = p_n80_y;
        projection_50_n80_z(i) = p_n80_z;
        
        projection_50_n120_x(i) = p_n120_x;
        projection_50_n120_y(i) = p_n120_y;
        projection_50_n120_z(i) = p_n120_z;

    end

    % VOI60
    for i = 1:length(original_x_60)
        [p_n0_x, p_n0_y, p_n0_z] = Orthproject(norm_0(1),norm_0(2),norm_0(3),...
            VOI60_centroid(1),VOI60_centroid(2),VOI60_centroid(3),...
            original_x_60(i),original_y_60(i),original_z_60(i));
        
        [p_n40_x, p_n40_y, p_n40_z] = Orthproject(norm_40(1),norm_40(2),norm_40(3),...
            VOI60_centroid(1),VOI60_centroid(2),VOI60_centroid(3),...
            original_x_60(i),original_y_60(i),original_z_60(i));
        
        [p_n80_x, p_n80_y, p_n80_z] = Orthproject(norm_80(1),norm_80(2),norm_80(3),...
            VOI60_centroid(1),VOI60_centroid(2),VOI60_centroid(3),...
            original_x_60(i),original_y_60(i),original_z_60(i));
        
        [p_n120_x, p_n120_y, p_n120_z] = Orthproject(norm_120(1),norm_120(2),norm_120(3),...
            VOI60_centroid(1),VOI60_centroid(2),VOI60_centroid(3),...
            original_x_60(i),original_y_60(i),original_z_60(i));
        
        projection_60_n0_x(i) = p_n0_x;
        projection_60_n0_y(i) = p_n0_y;
        projection_60_n0_z(i) = p_n0_z;
        
        projection_60_n40_x(i) = p_n40_x;
        projection_60_n40_y(i) = p_n40_y;
        projection_60_n40_z(i) = p_n40_z;
        
        projection_60_n80_x(i) = p_n80_x;
        projection_60_n80_y(i) = p_n80_y;
        projection_60_n80_z(i) = p_n80_z;
        
        projection_60_n120_x(i) = p_n120_x;
        projection_60_n120_y(i) = p_n120_y;
        projection_60_n120_z(i) = p_n120_z;
    end
    
    % VOI70
    for i = 1:length(original_x_70)
        [p_n0_x, p_n0_y, p_n0_z] = Orthproject(norm_0(1),norm_0(2),norm_0(3),...
            VOI70_centroid(1),VOI70_centroid(2),VOI70_centroid(3),...
            original_x_70(i),original_y_70(i),original_z_70(i));
        
        [p_n40_x, p_n40_y, p_n40_z] = Orthproject(norm_40(1),norm_40(2),norm_40(3),...
            VOI70_centroid(1),VOI70_centroid(2),VOI70_centroid(3),...
            original_x_70(i),original_y_70(i),original_z_70(i));
        
        [p_n80_x, p_n80_y, p_n80_z] = Orthproject(norm_80(1),norm_80(2),norm_80(3),...
            VOI70_centroid(1),VOI70_centroid(2),VOI70_centroid(3),...
            original_x_70(i),original_y_70(i),original_z_70(i));
        
        [p_n120_x, p_n120_y, p_n120_z] = Orthproject(norm_120(1),norm_120(2),norm_120(3),...
            VOI70_centroid(1),VOI70_centroid(2),VOI70_centroid(3),...
            original_x_70(i),original_y_70(i),original_z_70(i));
        
        projection_70_n0_x(i) = p_n0_x;
        projection_70_n0_y(i) = p_n0_y;
        projection_70_n0_z(i) = p_n0_z;
        
        projection_70_n40_x(i) = p_n40_x;
        projection_70_n40_y(i) = p_n40_y;
        projection_70_n40_z(i) = p_n40_z;
        
        projection_70_n80_x(i) = p_n80_x;
        projection_70_n80_y(i) = p_n80_y;
        projection_70_n80_z(i) = p_n80_z;
        
        projection_70_n120_x(i) = p_n120_x;
        projection_70_n120_y(i) = p_n120_y;
        projection_70_n120_z(i) = p_n120_z;
    end
        
    % VOI80
    for i = 1:length(original_x_80)
        [p_n0_x, p_n0_y, p_n0_z] = Orthproject(norm_0(1),norm_0(2),norm_0(3),...
            VOI80_centroid(1),VOI80_centroid(2),VOI80_centroid(3),...
            original_x_80(i),original_y_80(i),original_z_80(i));
        
        [p_n40_x, p_n40_y, p_n40_z] = Orthproject(norm_40(1),norm_40(2),norm_40(3),...
            VOI80_centroid(1),VOI80_centroid(2),VOI80_centroid(3),...
            original_x_80(i),original_y_80(i),original_z_80(i));
        
        [p_n80_x, p_n80_y, p_n80_z] = Orthproject(norm_80(1),norm_80(2),norm_80(3),...
            VOI80_centroid(1),VOI80_centroid(2),VOI80_centroid(3),...
            original_x_80(i),original_y_80(i),original_z_80(i));
        
        [p_n120_x, p_n120_y, p_n120_z] = Orthproject(norm_120(1),norm_120(2),norm_120(3),...
            VOI80_centroid(1),VOI80_centroid(2),VOI80_centroid(3),...
            original_x_80(i),original_y_80(i),original_z_80(i));
        
        projection_80_n0_x(i) = p_n0_x;
        projection_80_n0_y(i) = p_n0_y;
        projection_80_n0_z(i) = p_n0_z;
        
        projection_80_n40_x(i) = p_n40_x;
        projection_80_n40_y(i) = p_n40_y;
        projection_80_n40_z(i) = p_n40_z;
        
        projection_80_n80_x(i) = p_n80_x;
        projection_80_n80_y(i) = p_n80_y;
        projection_80_n80_z(i) = p_n80_z;
        
        projection_80_n120_x(i) = p_n120_x;
        projection_80_n120_y(i) = p_n120_y;
        projection_80_n120_z(i) = p_n120_z;
    end
    
    % VOI90
    for i = 1:length(original_x_90)
        
        [p_n0_x, p_n0_y, p_n0_z] = Orthproject(norm_0(1),norm_0(2),norm_0(3),...
            VOI90_centroid(1),VOI90_centroid(2),VOI90_centroid(3),...
            original_x_90(i),original_y_90(i),original_z_90(i));
        
        [p_n40_x, p_n40_y, p_n40_z] = Orthproject(norm_40(1),norm_40(2),norm_40(3),...
            VOI90_centroid(1),VOI90_centroid(2),VOI90_centroid(3),...
            original_x_90(i),original_y_90(i),original_z_90(i));
        
        [p_n80_x, p_n80_y, p_n80_z] = Orthproject(norm_80(1),norm_80(2),norm_80(3),...
            VOI90_centroid(1),VOI90_centroid(2),VOI90_centroid(3),...
            original_x_90(i),original_y_90(i),original_z_90(i));
        
        [p_n120_x, p_n120_y, p_n120_z] = Orthproject(norm_120(1),norm_120(2),norm_120(3),...
            VOI90_centroid(1),VOI90_centroid(2),VOI90_centroid(3),...
            original_x_90(i),original_y_90(i),original_z_90(i));
        
        projection_90_n0_x(i) = p_n0_x;
        projection_90_n0_y(i) = p_n0_y;
        projection_90_n0_z(i) = p_n0_z;
        
        projection_90_n40_x(i) = p_n40_x;
        projection_90_n40_y(i) = p_n40_y;
        projection_90_n40_z(i) = p_n40_z;
        
        projection_90_n80_x(i) = p_n80_x;
        projection_90_n80_y(i) = p_n80_y;
        projection_90_n80_z(i) = p_n80_z;
        
        projection_90_n120_x(i) = p_n120_x;
        projection_90_n120_y(i) = p_n120_y;
        projection_90_n120_z(i) = p_n120_z;
    end
    
    % VOI95
    for i = 1:length(original_x_95)
        
        [p_n0_x, p_n0_y, p_n0_z] = Orthproject(norm_0(1),norm_0(2),norm_0(3),...
            VOI95_centroid(1),VOI95_centroid(2),VOI95_centroid(3),...
            original_x_95(i),original_y_95(i),original_z_95(i));
        
        [p_n40_x, p_n40_y, p_n40_z] = Orthproject(norm_40(1),norm_40(2),norm_40(3),...
            VOI95_centroid(1),VOI95_centroid(2),VOI95_centroid(3),...
            original_x_95(i),original_y_95(i),original_z_95(i));
        
        [p_n80_x, p_n80_y, p_n80_z] = Orthproject(norm_80(1),norm_80(2),norm_80(3),...
            VOI95_centroid(1),VOI95_centroid(2),VOI95_centroid(3),...
            original_x_95(i),original_y_95(i),original_z_95(i));
        
        [p_n120_x, p_n120_y, p_n120_z] = Orthproject(norm_120(1),norm_120(2),norm_120(3),...
            VOI95_centroid(1),VOI95_centroid(2),VOI95_centroid(3),...
            original_x_95(i),original_y_95(i),original_z_95(i));
        
        projection_95_n0_x(i) = p_n0_x;
        projection_95_n0_y(i) = p_n0_y;
        projection_95_n0_z(i) = p_n0_z;
        
        projection_95_n40_x(i) = p_n40_x;
        projection_95_n40_y(i) = p_n40_y;
        projection_95_n40_z(i) = p_n40_z;
        
        projection_95_n80_x(i) = p_n80_x;
        projection_95_n80_y(i) = p_n80_y;
        projection_95_n80_z(i) = p_n80_z;
        
        projection_95_n120_x(i) = p_n120_x;
        projection_95_n120_y(i) = p_n120_y;
        projection_95_n120_z(i) = p_n120_z;

    end
    
    clearvars p_n0_x p_n0_y p_n0_z ...
        p_n40_x p_n40_y p_n40_z ...
        p_n80_x p_n80_y p_n80_z ...
        p_n120_x p_n120_y p_n120_z ...
        i
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change A to desired value
% TODO verander dit naar een functie
A = 50.40; % Notation A = VOI.Gantry_angle

switch A
    case 50.0       % VOI50, Gantry 0°
        % Perpendicular plane
        f1 = norm_0(1);
        f2 = norm_0(2);
        f3 = norm_0(3);
        f4 = -norm_0(1)*VOI50_centroid(1)-norm_0(2)*VOI50_centroid(2)-norm_0(3)*VOI50_centroid(3);
        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI50_centroid(1) + t * norm_0(1);
        Yl = VOI50_centroid(2) + t * norm_0(2);
        Zl = VOI50_centroid(3) + t * norm_0(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_50,original_y_50,original_z_50,'ob')
        % Projection
        scatter3(projection_50_n0_x,projection_50_n0_y,projection_50_n0_z,'or')
        hold off
        
    case 50.40      % VOI50, Gantry 40°
        % Perpendicular plane
        f1 = norm_40(1);
        f2 = norm_40(2);
        f3 = norm_40(3);
        f4 = -norm_40(1)*VOI50_centroid(1)-norm_40(2)*VOI50_centroid(2)-norm_40(3)*VOI50_centroid(3);

        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI50_centroid(1) + t * norm_40(1);
        Yl = VOI50_centroid(2) + t * norm_40(2);
        Zl = VOI50_centroid(3) + t * norm_40(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_50,original_y_50,original_z_50,'ob')
        % Projection
        scatter3(projection_50_n40_x,projection_50_n40_y,projection_50_n40_z,'or')
        hold off
        
    case 50.80      % VOI50, Gantry 80°
        % Perpendicular plane
        f1 = norm_80(1);
        f2 = norm_80(2);
        f3 = norm_80(3);
        f4 = -norm_80(1)*VOI50_centroid(1)-norm_80(2)*VOI50_centroid(2)-norm_80(3)*VOI50_centroid(3);

        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI50_centroid(1) + t * norm_80(1);
        Yl = VOI50_centroid(2) + t * norm_80(2);
        Zl = VOI50_centroid(3) + t * norm_80(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_50,original_y_50,original_z_50,'ob')
        % Projection
        scatter3(projection_50_n80_x,projection_50_n80_y,projection_50_n80_z,'or')
        hold off       
        
    case 50.120     % VOI50, Gantry 120°
        % Perpendicular plane
        f1 = norm_120(1);
        f2 = norm_120(2);
        f3 = norm_120(3);
        f4 = -norm_120(1)*VOI50_centroid(1)-norm_120(2)*VOI50_centroid(2)-norm_120(3)*VOI50_centroid(3);

        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI50_centroid(1) + t * norm_120(1);
        Yl = VOI50_centroid(2) + t * norm_120(2);
        Zl = VOI50_centroid(3) + t * norm_120(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_50,original_y_50,original_z_50,'ob')
        % Projection
        scatter3(projection_50_n120_x,projection_50_n120_y,projection_50_n120_z,'or')
        hold off 
        
    case 60.0       % VOI60, Gantry 0°
        % Perpendicular plane
        f1 = norm_0(1);
        f2 = norm_0(2);
        f3 = norm_0(3);
        f4 = -norm_0(1)*VOI60_centroid(1)-norm_0(2)*VOI60_centroid(2)-norm_0(3)*VOI60_centroid(3);
        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI60_centroid(1) + t * norm_0(1);
        Yl = VOI60_centroid(2) + t * norm_0(2);
        Zl = VOI60_centroid(3) + t * norm_0(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_60,original_y_60,original_z_60,'ob')
        % Projection
        scatter3(projection_60_n0_x,projection_60_n0_y,projection_60_n0_z,'or')
        hold off
        
    case 60.40      % VOI60, Gantry 40°
        % Perpendicular plane
        f1 = norm_40(1);
        f2 = norm_40(2);
        f3 = norm_40(3);
        f4 = -norm_40(1)*VOI60_centroid(1)-norm_40(2)*VOI60_centroid(2)-norm_40(3)*VOI60_centroid(3);
        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI60_centroid(1) + t * norm_40(1);
        Yl = VOI60_centroid(2) + t * norm_40(2);
        Zl = VOI60_centroid(3) + t * norm_40(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_60,original_y_60,original_z_60,'ob')
        % Projection
        scatter3(projection_60_n40_x,projection_60_n40_y,projection_60_n40_z,'or')
        hold off
        
    case 60.80      % VOI60, Gantry 80°
        % Perpendicular plane
        f1 = norm_80(1);
        f2 = norm_80(2);
        f3 = norm_80(3);
        f4 = -norm_80(1)*VOI60_centroid(1)-norm_80(2)*VOI60_centroid(2)-norm_80(3)*VOI60_centroid(3);
        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI60_centroid(1) + t * norm_80(1);
        Yl = VOI60_centroid(2) + t * norm_80(2);
        Zl = VOI60_centroid(3) + t * norm_80(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_60,original_y_60,original_z_60,'ob')
        % Projection
        scatter3(projection_60_n80_x,projection_60_n80_y,projection_60_n80_z,'or')
        hold off
        
    case 60.120     % VOI60, Gantry 120°
        % Perpendicular plane
        f1 = norm_120(1);
        f2 = norm_120(2);
        f3 = norm_120(3);
        f4 = -norm_120(1)*VOI60_centroid(1)-norm_120(2)*VOI60_centroid(2)-norm_120(3)*VOI60_centroid(3);
        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI60_centroid(1) + t * norm_120(1);
        Yl = VOI60_centroid(2) + t * norm_120(2);
        Zl = VOI60_centroid(3) + t * norm_120(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_60,original_y_60,original_z_60,'ob')
        % Projection
        scatter3(projection_60_n120_x,projection_60_n120_y,projection_60_n120_z,'or')
        hold off
        
    case 70.0       % VOI70, Gantry 0°
        % Perpendicular plane
        f1 = norm_0(1);
        f2 = norm_0(2);
        f3 = norm_0(3);
        f4 = -norm_0(1)*VOI70_centroid(1)-norm_0(2)*VOI70_centroid(2)-norm_0(3)*VOI70_centroid(3);
        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI70_centroid(1) + t * norm_0(1);
        Yl = VOI70_centroid(2) + t * norm_0(2);
        Zl = VOI70_centroid(3) + t * norm_0(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_70,original_y_70,original_z_70,'ob')
        % Projection
        scatter3(projection_70_n0_x,projection_70_n0_y,projection_70_n0_z,'or')
        hold off
        
    case 70.40      % VOI70, Gantry 40°
        % Perpendicular plane
        f1 = norm_40(1);
        f2 = norm_40(2);
        f3 = norm_40(3);
        f4 = -norm_40(1)*VOI70_centroid(1)-norm_40(2)*VOI70_centroid(2)-norm_40(3)*VOI70_centroid(3);
        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI70_centroid(1) + t * norm_40(1);
        Yl = VOI70_centroid(2) + t * norm_40(2);
        Zl = VOI70_centroid(3) + t * norm_40(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_70,original_y_70,original_z_70,'ob')
        % Projection
        scatter3(projection_70_n40_x,projection_70_n40_y,projection_70_n40_z,'or')
        hold off
        
    case 70.80      % VOI70, Gantry 80°
        % Perpendicular plane
        f1 = norm_80(1);
        f2 = norm_80(2);
        f3 = norm_80(3);
        f4 = -norm_80(1)*VOI70_centroid(1)-norm_80(2)*VOI70_centroid(2)-norm_80(3)*VOI70_centroid(3);
        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI70_centroid(1) + t * norm_80(1);
        Yl = VOI70_centroid(2) + t * norm_80(2);
        Zl = VOI70_centroid(3) + t * norm_80(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_70,original_y_70,original_z_70,'ob')
        % Projection
        scatter3(projection_70_n80_x,projection_70_n80_y,projection_70_n80_z,'or')
        hold off
        
    case 70.120     % VOI70, Gantry 120°
        % Perpendicular plane
        f1 = norm_120(1);
        f2 = norm_120(2);
        f3 = norm_120(3);
        f4 = -norm_120(1)*VOI70_centroid(1)-norm_120(2)*VOI70_centroid(2)-norm_120(3)*VOI70_centroid(3);
        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI70_centroid(1) + t * norm_120(1);
        Yl = VOI70_centroid(2) + t * norm_120(2);
        Zl = VOI70_centroid(3) + t * norm_120(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_70,original_y_70,original_z_70,'ob')
        % Projection
        scatter3(projection_70_n120_x,projection_70_n120_y,projection_70_n120_z,'or')
        hold off
        
    case 80.0
    case 80.40
    case 80.80
    case 80.120
    case 90.0
    case 90.40
    case 90.80
    case 90.120
    case 95.0
        
        % Perpendicular plane
        f1 = norm_0(1);
        f2 = norm_0(2);
        f3 = norm_0(3);
        f4 = -norm_0(1)*VOI95_centroid(1)-norm_0(2)*VOI95_centroid(2)-norm_0(3)*VOI95_centroid(3);
        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI95_centroid(1) + t * norm_0(1);
        Yl = VOI95_centroid(2) + t * norm_0(2);
        Zl = VOI95_centroid(3) + t * norm_0(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_95,original_y_95,original_z_95,'ob')
        % Projection
        scatter3(projection_95_n0_x,projection_95_n0_y,projection_95_n0_z,'or')
        hold off
        
    case 95.40
        % Perpendicular plane
        f1 = norm_40(1);
        f2 = norm_40(2);
        f3 = norm_40(3);
        f4 = -norm_40(1)*VOI95_centroid(1)-norm_40(2)*VOI95_centroid(2)-norm_40(3)*VOI95_centroid(3);
        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI95_centroid(1) + t * norm_40(1);
        Yl = VOI95_centroid(2) + t * norm_40(2);
        Zl = VOI95_centroid(3) + t * norm_40(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_95,original_y_95,original_z_95,'ob')
        % Projection
        scatter3(projection_95_n40_x,projection_95_n40_y,projection_95_n40_z,'or')
        hold off
        
    case 95.80
        % Perpendicular plane
        f1 = norm_80(1);
        f2 = norm_80(2);
        f3 = norm_80(3);
        f4 = -norm_80(1)*VOI95_centroid(1)-norm_80(2)*VOI95_centroid(2)-norm_80(3)*VOI95_centroid(3);
        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI95_centroid(1) + t * norm_80(1);
        Yl = VOI95_centroid(2) + t * norm_80(2);
        Zl = VOI95_centroid(3) + t * norm_80(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_95,original_y_95,original_z_95,'ob')
        % Projection
        scatter3(projection_95_n80_x,projection_95_n80_y,projection_95_n80_z,'or')
        hold off
        
    case 95.120
        % Perpendicular plane
        f1 = norm_120(1);
        f2 = norm_120(2);
        f3 = norm_120(3);
        f4 = -norm_120(1)*VOI95_centroid(1)-norm_120(2)*VOI95_centroid(2)-norm_120(3)*VOI95_centroid(3);
        [x, y] = meshgrid(1:5:500); % Generate x and y data
        z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
        surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
        hold on
        xlabel('X-as')
        ylabel('Y-as')
        zlabel('Z-as')
        xlim([0 500]) 
        ylim([0 500])
        zlim([0 500])

        % Incident beam
        t = -1000:1:1000;
        Xl = VOI95_centroid(1) + t * norm_120(1);
        Yl = VOI95_centroid(2) + t * norm_120(2);
        Zl = VOI95_centroid(3) + t * norm_120(3);
        plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

        % Tumour plot
        % Tumour
        scatter3(original_x_95,original_y_95,original_z_95,'ob')
        % Projection
        scatter3(projection_95_n120_x,projection_95_n120_y,projection_95_n120_z,'or')
        hold off
        
    otherwise
        
        warning('Not a valid option')
        
end

clearvars A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BOUNDING BOX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Declaration
    t_VOI50 = zeros(16,3);
    t_VOI60 = zeros(16,3);
    t_VOI70 = zeros(16,3);
    t_VOI80 = zeros(16,3);
    t_VOI90 = zeros(16,3);
    t_VOI95 = zeros(16,3);

% Gantry 0°
    f1 = norm_0(1);
    f2 = norm_0(2);
    f3 = norm_0(3);

    % VOI50
    f4 = -norm_0(1)*VOI50_centroid(1)-norm_0(2)*VOI50_centroid(2)-norm_0(3)*VOI50_centroid(3);
    pro_min_y = min(projection_50_n0_y);
    pro_max_y = max(projection_50_n0_y);
    pro_min_z = min(projection_50_n0_z);
    pro_max_z = max(projection_50_n0_z);  

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI50(1:4,1:3) = Pts;

    % VOI60
    f4 = -norm_0(1)*VOI60_centroid(1)-norm_0(2)*VOI60_centroid(2)-norm_0(3)*VOI60_centroid(3);
    pro_min_y = min(projection_60_n0_y);
    pro_max_y = max(projection_60_n0_y);
    pro_min_z = min(projection_60_n0_z);
    pro_max_z = max(projection_60_n0_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI60(1:4,1:3) = Pts;
    
    % VOI70
    f4 = -norm_0(1)*VOI70_centroid(1)-norm_0(2)*VOI70_centroid(2)-norm_0(3)*VOI70_centroid(3);
    pro_min_y = min(projection_70_n0_y);
    pro_max_y = max(projection_70_n0_y);
    pro_min_z = min(projection_70_n0_z);
    pro_max_z = max(projection_70_n0_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI70(1:4,1:3) = Pts;
    
    % VOI80
    f4 = -norm_0(1)*VOI80_centroid(1)-norm_0(2)*VOI80_centroid(2)-norm_0(3)*VOI80_centroid(3);
    pro_min_y = min(projection_80_n0_y);
    pro_max_y = max(projection_80_n0_y);
    pro_min_z = min(projection_80_n0_z);
    pro_max_z = max(projection_80_n0_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI80(1:4,1:3) = Pts;
    
    % VOI90
    f4 = -norm_0(1)*VOI90_centroid(1)-norm_0(2)*VOI90_centroid(2)-norm_0(3)*VOI90_centroid(3);
    pro_min_y = min(projection_90_n0_y);
    pro_max_y = max(projection_90_n0_y);
    pro_min_z = min(projection_90_n0_z);
    pro_max_z = max(projection_90_n0_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI90(1:4,1:3) = Pts;
    
    % VOI95
    f4 = -norm_0(1)*VOI95_centroid(1)-norm_0(2)*VOI95_centroid(2)-norm_0(3)*VOI95_centroid(3);
    pro_min_y = min(projection_95_n0_y);
    pro_max_y = max(projection_95_n0_y);
    pro_min_z = min(projection_95_n0_z);
    pro_max_z = max(projection_95_n0_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI95(1:4,1:3) = Pts;
    
% Gantry 40°
    f1 = norm_40(1);
    f2 = norm_40(2);
    f3 = norm_40(3);

    % VOI50
    f4 = -norm_40(1)*VOI50_centroid(1)-norm_40(2)*VOI50_centroid(2)-norm_40(3)*VOI50_centroid(3);
    pro_min_y = min(projection_50_n40_y);
    pro_max_y = max(projection_50_n40_y);
    pro_min_z = min(projection_50_n40_z);
    pro_max_z = max(projection_50_n40_z);
    
    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI50(5:8,1:3) = Pts;

    % VOI60
    f4 = -norm_40(1)*VOI60_centroid(1)-norm_40(2)*VOI60_centroid(2)-norm_40(3)*VOI60_centroid(3);
    pro_min_y = min(projection_60_n40_y);
    pro_max_y = max(projection_60_n40_y);
    pro_min_z = min(projection_60_n40_z);
    pro_max_z = max(projection_60_n40_z); 
    
    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI60(5:8,1:3) = Pts;

    % VOI70
    f4 = -norm_40(1)*VOI70_centroid(1)-norm_40(2)*VOI70_centroid(2)-norm_40(3)*VOI70_centroid(3);
    pro_min_y = min(projection_70_n40_y);
    pro_max_y = max(projection_70_n40_y);
    pro_min_z = min(projection_70_n40_z);
    pro_max_z = max(projection_70_n40_z);
    
    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI70(5:8,1:3) = Pts;

    % VOI80
    f4 = -norm_40(1)*VOI80_centroid(1)-norm_40(2)*VOI80_centroid(2)-norm_40(3)*VOI80_centroid(3);
    pro_min_y = min(projection_80_n40_y);
    pro_max_y = max(projection_80_n40_y);
    pro_min_z = min(projection_80_n40_z);
    pro_max_z = max(projection_80_n40_z);
    
    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI80(5:8,1:3) = Pts;

    % VOI90
    f4 = -norm_40(1)*VOI90_centroid(1)-norm_40(2)*VOI90_centroid(2)-norm_40(3)*VOI90_centroid(3);
    pro_min_y = min(projection_90_n40_y);
    pro_max_y = max(projection_90_n40_y);
    pro_min_z = min(projection_90_n40_z);
    pro_max_z = max(projection_90_n40_z);
    
    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI90(5:8,1:3) = Pts;

    % VOI95
    f4 = -norm_40(1)*VOI95_centroid(1)-norm_40(2)*VOI95_centroid(2)-norm_40(3)*VOI95_centroid(3);
    pro_min_y = min(projection_95_n40_y);
    pro_max_y = max(projection_95_n40_y);
    pro_min_z = min(projection_95_n40_z);
    pro_max_z = max(projection_95_n40_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI95(5:8,1:3) = Pts;
    
% Gantry 80°
    f1 = norm_80(1);
    f2 = norm_80(2);
    f3 = norm_80(3);

    % VOI50
    f4 = -norm_80(1)*VOI50_centroid(1)-norm_80(2)*VOI50_centroid(2)-norm_80(3)*VOI50_centroid(3);
    pro_min_y = min(projection_50_n80_y);
    pro_max_y = max(projection_50_n80_y);
    pro_min_z = min(projection_50_n80_z);
    pro_max_z = max(projection_50_n80_z);
    
    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI50(9:12,1:3) = Pts;

    % VOI60
    f4 = -norm_80(1)*VOI60_centroid(1)-norm_80(2)*VOI60_centroid(2)-norm_80(3)*VOI60_centroid(3);
    pro_min_y = min(projection_60_n80_y);
    pro_max_y = max(projection_60_n80_y);
    pro_min_z = min(projection_60_n80_z);
    pro_max_z = max(projection_60_n80_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI60(9:12,1:3) = Pts;
    
    % VOI70
    f4 = -norm_80(1)*VOI70_centroid(1)-norm_80(2)*VOI70_centroid(2)-norm_80(3)*VOI70_centroid(3);
    pro_min_y = min(projection_70_n80_y);
    pro_max_y = max(projection_70_n80_y);
    pro_min_z = min(projection_70_n80_z);
    pro_max_z = max(projection_70_n80_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI70(9:12,1:3) = Pts;
    
    % VOI80
    f4 = -norm_80(1)*VOI80_centroid(1)-norm_80(2)*VOI80_centroid(2)-norm_80(3)*VOI80_centroid(3);
    pro_min_y = min(projection_80_n80_y);
    pro_max_y = max(projection_80_n80_y);
    pro_min_z = min(projection_80_n80_z);
    pro_max_z = max(projection_80_n80_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI80(9:12,1:3) = Pts;    
    
    % VOI90
    f4 = -norm_80(1)*VOI90_centroid(1)-norm_80(2)*VOI90_centroid(2)-norm_80(3)*VOI90_centroid(3);
    pro_min_y = min(projection_90_n80_y);
    pro_max_y = max(projection_90_n80_y);
    pro_min_z = min(projection_90_n80_z);
    pro_max_z = max(projection_90_n80_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI90(9:12,1:3) = Pts;
    
    % VOI95
    f4 = -norm_80(1)*VOI95_centroid(1)-norm_80(2)*VOI95_centroid(2)-norm_80(3)*VOI95_centroid(3);
    pro_min_y = min(projection_95_n80_y);
    pro_max_y = max(projection_95_n80_y);
    pro_min_z = min(projection_95_n80_z);
    pro_max_z = max(projection_95_n80_z);
    
    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI95(9:12,1:3) = Pts;

% Gantry 120°
    f1 = norm_120(1);
    f2 = norm_120(2);
    f3 = norm_120(3);

    % VOI50
    f4 = -norm_120(1)*VOI50_centroid(1)-norm_120(2)*VOI50_centroid(2)-norm_120(3)*VOI50_centroid(3);
    pro_min_y = min(projection_50_n120_y);
    pro_max_y = max(projection_50_n120_y);
    pro_min_z = min(projection_50_n120_z);
    pro_max_z = max(projection_50_n120_z);  

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI50(13:16,1:3) = Pts;
    
    % VOI60
    f4 = -norm_120(1)*VOI60_centroid(1)-norm_120(2)*VOI60_centroid(2)-norm_120(3)*VOI60_centroid(3);
    pro_min_y = min(projection_60_n120_y);
    pro_max_y = max(projection_60_n120_y);
    pro_min_z = min(projection_60_n120_z);
    pro_max_z = max(projection_60_n120_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI60(13:16,1:3) = Pts;
    
    % VOI70
    f4 = -norm_120(1)*VOI70_centroid(1)-norm_120(2)*VOI70_centroid(2)-norm_120(3)*VOI70_centroid(3);
    pro_min_y = min(projection_70_n120_y);
    pro_max_y = max(projection_70_n120_y);
    pro_min_z = min(projection_70_n120_z);
    pro_max_z = max(projection_70_n120_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI70(13:16,1:3) = Pts;
    
    % VOI80
    f4 = -norm_120(1)*VOI80_centroid(1)-norm_120(2)*VOI80_centroid(2)-norm_120(3)*VOI80_centroid(3);
    pro_min_y = min(projection_80_n120_y);
    pro_max_y = max(projection_80_n120_y);
    pro_min_z = min(projection_80_n120_z);
    pro_max_z = max(projection_80_n120_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI80(13:16,1:3) = Pts;
    
    % VOI90
    f4 = -norm_120(1)*VOI90_centroid(1)-norm_120(2)*VOI90_centroid(2)-norm_120(3)*VOI90_centroid(3);
    pro_min_y = min(projection_90_n120_y);
    pro_max_y = max(projection_90_n120_y);
    pro_min_z = min(projection_90_n120_z);
    pro_max_z = max(projection_90_n120_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI90(13:16,1:3) = Pts;
    
    % VOI95
    f4 = -norm_120(1)*VOI95_centroid(1)-norm_120(2)*VOI95_centroid(2)-norm_120(3)*VOI95_centroid(3);
    pro_min_y = min(projection_95_n120_y);
    pro_max_y = max(projection_95_n120_y);
    pro_min_z = min(projection_95_n120_z);
    pro_max_z = max(projection_95_n120_z); 

    Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
    X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
    Pts(:,1) = X;
    t_VOI95(13:16,1:3) = Pts;
    
% Generate table
    Gantry = [0; 0; 0; 0; 40; 40; 40; 40; 80; 80; 80; 80; 120; 120; 120; 120];
    BoundingBoxPoints = table(Gantry,t_VOI50,t_VOI60,t_VOI70,Gantry,t_VOI80,t_VOI90,t_VOI95,Gantry);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT BOUNDING BOX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change A to desired value
% TODO verander dit naar een functie
A = 50.0; % Notation A = VOI.Gantry_angle

% switch A
%     case 50.0       % VOI50, Gantry 0°
%     case 50.40      % VOI50, Gantry 40°
%     case 50.80      % VOI50, Gantry 80°
%     case 50.120     % VOI50, Gantry 120°
%     case 60.0       % VOI60, Gantry 0°
%     case 60.40      % VOI60, Gantry 40°
%     case 60.80      % VOI60, Gantry 80°
%     case 60.120     % VOI60, Gantry 120°
%     case 70.0       % VOI70, Gantry 0°
%     case 70.40      % VOI70, Gantry 40°
%     case 70.80      % VOI70, Gantry 80°
%     case 70.120     % VOI70, Gantry 120°
%     case 80.0
%     case
%     case
%     case
%     otherwise
%         
%         error('This is not a valid input')
%         
% end

clearvars A




figure()
f1 = norm_80(1);
f2 = norm_80(2);
f3 = norm_80(3);
f4 = -norm_80(1)*VOI50_centroid(1)-norm_80(2)*VOI50_centroid(2)-norm_80(3)*VOI50_centroid(3);

[x, y] = meshgrid(1:5:500); % Generate x and y data
z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
surf(x,y,z,'EdgeColor','k','FaceColor','none') %Plot the surface
hold on
xlabel('X-as')
ylabel('Y-as')
zlabel('Z-as')
xlim([0 500]) 
ylim([0 500])
zlim([0 500])

% Incident beam
t = -1000:1:1000;
Xl = VOI50_centroid(1) + t * norm_80(1);
Yl = VOI50_centroid(2) + t * norm_80(2);
Zl = VOI50_centroid(3) + t * norm_80(3);
plot3(Xl,Yl,Zl,'.k') % Line through isocenter, orthogonal to projection plane

% Tumour plot
% Tumour
scatter3(original_x_50,original_y_50,original_z_50,'ob')
% Projection
scatter3(projection_50_n80_x,projection_50_n80_y,projection_50_n80_z,'or')
        
pro_min_x = min(projection_50_n80_x);
pro_max_x = max(projection_50_n80_x);
pro_min_y = min(projection_50_n80_y);
pro_max_y = max(projection_50_n80_y);
pro_min_z = min(projection_50_n80_z);
pro_max_z = max(projection_50_n80_z);  

images.roi.Cuboid(gca,'Color','black','InteractionsAllowed','none',...
    'LineWidth',2,'Position',[pro_min_x pro_min_y pro_min_z ...
    pro_max_x-pro_min_x pro_max_y-pro_min_y pro_max_z-pro_min_z]);


% Perpendicular plane
f1 = norm_80(1);
f2 = norm_80(2);
f3 = norm_80(3);
f4 = -norm_80(1)*VOI50_centroid(1)-norm_80(2)*VOI50_centroid(2)-norm_80(3)*VOI50_centroid(3);

Pts = [0 pro_min_y pro_min_z; 0 pro_min_y pro_max_z; 0 pro_max_y pro_min_z; 0 pro_max_y pro_max_z];
% P1 Ymin Zmin
% P2 Ymin Zmax
% P3 Ymax Zmin
% P4 Ymax Zmax TODO telkens Y en Z gebruiken??

X = -(f2.*Pts(:,2)+f3.*Pts(:,3)+f4)./f1;
Pts(:,1) = X;

scatter3(gca,Pts(:,1),Pts(:,2),Pts(:,3),'cyan')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% The following shows two methods for constructing, an arbitrary plane of the form:
%  Ax + By + Cz + D = 0,
% where the coefficients "A", "B", "C", and "D" are known values.
% Method 1: The PATCH Function


% (1) Equation plane
% a*(x-x1) + b*(y-y1) + c*(z-z1) = 0
% (2) Equation Plane 
% f1*x + f2*y + f3*z + f4 = 0
% f1 = a
% f2 = b
% f3 = c
% f4 = -a*x1-b*y1-c*z1

f1 = norm_0(1);
f2 = norm_0(2);
f3 = norm_0(3);
f4 = -norm_0(1)*VOI50_centroid(1)-norm_0(2)*VOI50_centroid(2)-norm_0(3)*VOI50_centroid(3);

% plot plane
x = [1000 -1000 -1000 1000]; % Generate data for x vertices
y = [1000 1000 -1000 -1000]; % Generate data for y vertices
z = zeros(500,500); %-1/f3*(f1*x + f2*y + f4); % Solve for z vertices data
patch(x, y, z);
view(3)

xlabel('X-as')
ylabel('Y-as')
zlabel('Z-as')
xlim([0 500]) 
ylim([0 500])
zlim([0 500])


figure()
[x, y] = meshgrid(1:1:500); % Generate x and y data
z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
%z = zeros(500,500);
surf(x,y,z,'FaceColor','k') %Plot the surface
xlabel('X-as')
ylabel('Y-as')
zlabel('Z-as')
xlim([0 500]) 
ylim([0 500])
zlim([0 500])




Point_ex = VOI50_centroid;
Vector_ex = [0, 0, 1];
w = null(Vector_ex); % Find two orthonormal vectors which are orthogonal to v
[P,Q] = meshgrid(0:5:500); % Provide a gridwork (you choose the size)
X = Point_ex(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
Y = Point_ex(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
Z = Point_ex(3)+w(3,1)*P+w(3,2)*Q;
surf(X,Y,Z,'EdgeColor','k','FaceColor','k')
xlabel('X-as')
ylabel('Y-as')
zlabel('Z-as')
xlim([0 500]) 
ylim([0 500])
zlim([0 500])
hold on
plot3(original(:,1),original(:,2),original(:,3),'ob')
plot3(projection(:,1),projection(:,2),projection(:,3),'or')

Vector_ex = [1,0,1/tand(40)];
w = null(Vector_ex); % Find two orthonormal vectors which are orthogonal to v
[P,Q] = meshgrid(0:5:500); % Provide a gridwork (you choose the size)
X = Point_ex(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
Y = Point_ex(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
Z = Point_ex(3)+w(3,1)*P+w(3,2)*Q;
surf(X,Y,Z,'EdgeColor','b','FaceColor','b')

Vector_ex = [1,0,1/tand(80)];
w = null(Vector_ex); % Find two orthonormal vectors which are orthogonal to v
[P,Q] = meshgrid(0:5:500); % Provide a gridwork (you choose the size)
X = Point_ex(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
Y = Point_ex(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
Z = Point_ex(3)+w(3,1)*P+w(3,2)*Q;
surf(X,Y,Z,'EdgeColor','g','FaceColor','g')

Vector_ex = [1,0,tand(30)];
w = null(Vector_ex); % Find two orthonormal vectors which are orthogonal to v
[P,Q] = meshgrid(0:5:500); % Provide a gridwork (you choose the size)
X = Point_ex(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
Y = Point_ex(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
Z = Point_ex(3)+w(3,1)*P+w(3,2)*Q;
surf(X,Y,Z,'EdgeColor','m','FaceColor','m')
hold off



% Extract image information
% NOTE 
% It's [left, top, width, height]. 
% Just be aware that the left and top are 0.5 pixels to the left and above, respectively, 
% than the actual first column and first row of the binary image. 
% It's because the Bounding box is defined such that it "contains" the blob, 
% and if it ran right through the center of the top-most and left-most pixel, 
% then that definition becomes somewhat ambiguous.


% I = imread('someColorImage.png'); %load an image
% BW = I(:,:,1)+I(:,:,2)+I(:,:,3); %make it into a black/white img
% s = regionprops(BW,'BoundingBox'); %get individual boundingboxes
% s = cat(1, s.BoundingBox); %concatenate struct into matrix
% %s has structure [left bottom width height]
% t = [s(:,1) s(:,2) s(:,1)+s(:,3) s(:,2)+s(:,4)]; 
% %t has structure [left bottom right top]
% leftBottomCorner = min(t(:,1:2));
% rightTopCorner = max(t(:,3:4));
% %boundingbox, again in [left bottom width height] style
% bb = [leftBottomCorner rightTopCorner-leftBottomCorner];

trisurf(box,original(:,1),original(:,2),original(:,3),'Facecolor','red','FaceAlpha',0.1)
trisurf(box2,projection(:,1),projection(:,2),projection(:,3),'Facecolor','red','FaceAlpha',0.1)

% (1)
     x = rand(100,1);
     y = rand(100,1);
     z = rand(100,1);
     [rotmat,cornerpoints,volume,surface] = minboundbox(x,y,z);
     plot3(x,y,z,'b.');hold on;plotminbox(cornerpoints,'r');
     
% (2)
     x=[1,0,0.1,1];y=[0,1,0.1,1];z=[0,0,0.9,1];
     [nerd,cornerpoints1,nerd,nerd] = minboundbox(x,y,z,'v',1);
     [nerd,cornerpoints2,nerd,nerd] = minboundbox(x,y,z,'v',2);
     [nerd,cornerpoints3,nerd,nerd] = minboundbox(x,y,z,'v',3);
     plot3(x,y,z,'bo','LineWidth',5);hold on;
     plotminbox(cornerpoints1,'b');hold on;
     plotminbox(cornerpoints2,'m');hold on;
     plotminbox(cornerpoints3,'r');axis equal;grid on;hold off;

% example

    
    [~,cornerpoints,~,~] = minboundbox(projection(:,1),projection(:,2),projection(:,3),'v',4)
    plot3(original(:,1),original(:,2),original(:,3),'or')
    hold on
    plot3(projection(:,1),projection(:,2),projection(:,3),'og')
    plotminbox(cornerpoints,'b')
    
     
     
box2 = boundary(projection(:,1),projection(:,2),projection(:,3),1)


% multiple bounding boxes to one

% Convert the [x y width height] coordinates to start and end coordinates,
% [x y x+width-1 y+height-1]
% Now take the minimum of all of the left x over all of the boxes to get the left bound, the maximum over all of the right x to get the right bound, the minimum over the top y to get the lower bound, the maximum over the bottom y to get the upper bound.
% 

BoundingBox = images.roi.Cuboid(example,'Color','black','InteractionsAllowed','none',...
    'LineWidth',2,'Position',[pro_min_x pro_min_y pro_min_z pro_max_x-pro_min_x pro_max_y-pro_min-y pro_max_z-pro_min_z]);


% figure()
% g = imshow(PET_VOI_SEG_LARG_50(:,:,85))
% 
% figure()
% h = imagesc(PET_VOI_SEG_LARG_50(:,:,85))
% [x y] = getpts

% for i = 1:size(image,1)
%     for j = 1:size(image,2)
%         for k = 1:size(image,3)
%             
%             value = PET_VOI_SEG_LARG_50(i,j,k);
%             
%             if value
%                 [px, py, pz] = Orthproject(norm_0(1),norm_0(2),norm_0(3),VOI50_centroid(1),VOI50_centroid(2),VOI50_centroid(3),i,j,k);
%                 
%                 original(counter,1) = i; % omgewisseld
%                 original(counter,2) = j; % omgewisseld
%                 original(counter,3) = k;
%                 
%                 projection(counter,1) = px; % omgewisseld
%                 projection(counter,2) = py; % omgewisseld
%                 projection(counter,3) = pz;
%                 
%                 counter = counter + 1;
%             end
%         end
%     end
% end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT ALL PLANES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()

% Gantry 0°
    f1 = norm_0(1);
    f2 = norm_0(2);
    f3 = norm_0(3);
    f4 = -norm_0(1)*VOI50_centroid(1)-norm_0(2)*VOI50_centroid(2)-norm_0(3)*VOI50_centroid(3);

    [x, y] = meshgrid(1:1:500); % Generate x and y data
    z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
    surf(x,y,z,'EdgeColor','k','FaceColor','k') %Plot the surface
    hold on
    xlabel('X-as')
    ylabel('Y-as')
    zlabel('Z-as')
    xlim([0 500]) 
    ylim([0 500])
    zlim([0 500])

    % Incident beam
        t = -1000:1:1000;
        Xl = VOI50_centroid(1) + t * norm_0(1);
        Yl = VOI50_centroid(2) + t * norm_0(2);
        Zl = VOI50_centroid(3) + t * norm_0(3);
        l1 = plot3(Xl,Yl,Zl,'.k')                                                % line through isocenter
        
% Gantry 40°
    f1 = norm_40(1);
    f2 = norm_40(2);
    f3 = norm_40(3);
    f4 = -norm_40(1)*VOI50_centroid(1)-norm_40(2)*VOI50_centroid(2)-norm_40(3)*VOI50_centroid(3);

    [x, y] = meshgrid(1:1:500); % Generate x and y data
    z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
    surf(x,y,z,'EdgeColor','r','FaceColor','r') %Plot the surface
    
    % Incident beam
        t = -1000:1:1000;
        Xl = VOI50_centroid(1) + t * norm_40(1);
        Yl = VOI50_centroid(2) + t * norm_40(2);
        Zl = VOI50_centroid(3) + t * norm_40(3);
        l2 = plot3(Xl,Yl,Zl,'.r')                                                % line through isocenter
        

% Gantry 80°
    f1 = norm_80(1);
    f2 = norm_80(2);
    f3 = norm_80(3);
    f4 = -norm_80(1)*VOI50_centroid(1)-norm_80(2)*VOI50_centroid(2)-norm_80(3)*VOI50_centroid(3);

    [x, y] = meshgrid(1:1:500); % Generate x and y data
    z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
    surf(x,y,z,'EdgeColor','b','FaceColor','b') %Plot the surface

    % Incident beam
        t = -1000:1:1000;
        Xl = VOI50_centroid(1) + t * norm_80(1);
        Yl = VOI50_centroid(2) + t * norm_80(2);
        Zl = VOI50_centroid(3) + t * norm_80(3);
        l3 = plot3(Xl,Yl,Zl,'.b')                                                % line through isocenter


% Gantry 120°
    f1 = norm_120(1);
    f2 = norm_120(2);
    f3 = norm_120(3);
    f4 = -norm_120(1)*VOI50_centroid(1)-norm_120(2)*VOI50_centroid(2)-norm_120(3)*VOI50_centroid(3);

    [x, y] = meshgrid(1:1:500); % Generate x and y data
    z = -1./f3*(f1.*x + f2.*y + f4); % Solve for z data
    surf(x,y,z,'EdgeColor','g','FaceColor','g') %Plot the surface
    
    % Incident beam
        t = -1000:1:1000;
        Xl = VOI50_centroid(1) + t * norm_120(1);
        Yl = VOI50_centroid(2) + t * norm_120(2);
        Zl = VOI50_centroid(3) + t * norm_120(3);
        l4 = plot3(Xl,Yl,Zl,'.g')                                                % line through isocenter

% volume
scatter3(original_x_50,original_y_50,original_z_50,'ob')                    % original tumour volume
scatter3(projection_50_n0_x,projection_50_n0_y,projection_50_n0_z,'or')              % projected tumour volume
scatter3(VOI50_centroid(1),VOI50_centroid(2),VOI50_centroid(3),'*y')
legend([l1, l2, l3, l4],{'Gantry 0°','Gantry 40°', 'Gantry 80°', 'Gantry 120°'})
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% JSON ENCODING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

json.FileVersion = 5;
json.name = "NewExperiment_SubjectHierarchy";
json.studies.beams.adjust_weight = "Yes";
json.studies.beams.collimator = "v5.0x5.0 mm";
json.studies.beams.collimator_rotation_angle = "0";
json.studies.beams.collimator_x = 5.0; 
json.studies.beams.collimator_y = 5.0; 
json.studies.beams.couch = 0; 
json.studies.beams.gantry = 0; 
json.studies.beams.isocenter = "IsoC_1"; 
json.studies.beams.label = "Beam"; 
json.studies.beams.ssd = "null"; %eigelijk zonder aanhalingstekens
json.studies.beams.time = "0"; 
json.studies.beams.type = "Beam"; 
json.studies.beams.weight = 25.0;
json.studies.contours = [];
json.studies.dose_type = "Water";
json.studies.heterogene = "true";
json.studies.isocenters.coordinates = [0.5, 11.4, -6.26];
json.studies.isocenters.dose = 2000;
json.studies.isocenters.name = "IsoC_1";

encoded = jsonencode(json)

FileID = fopen('Test.json','w')

fprintf(FileID,encoded)
