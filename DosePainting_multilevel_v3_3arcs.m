%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Dose Painting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%RADIATION PLANNING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information
%
% This script calculates the setup parameters for PET based radiation 
% therapy with the SARRP (XStrahl). 
%
% Specifications/limitations SARRP:
%       - Couch positions:  0°, -45°, -90°
%       - Gantry angles:    any angle
%                           for 0° & -45°   -> 0° - 120°
%                           for -90°        -> 0° - 60° (avoid collision
%                           with animal)
%       - Jaw size:         maximum: 40 x 80 mm 
%                           minmum: 1 x 1 mm
%       - SARRP resolution  yaw size = 0.01 mm
%                           isocenter = 0.01 mm
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
disp('Reading Toolboxes...')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Preclinical\Scripts\App_BoundingBox')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Preclinical\DosePainting\')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VARIABLES TO ADJUST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% beams
angles      = [0, 20, 40, 60, 80, 100, 120, 0, 20, 40, 60, 80, 100, 120, 0, 20, 40, 60];    % Gantry angles
couch_pos   = [0, 0, 0, 0, 0, 0, 0, -45, -45, -45, -45, -45, -45, -45, -90, -90, -90, -90]; % Couch positions

if length(angles) ~= length(couch_pos)
    disp('Variables ''angles'' and ''couch_pos'' should have the same length.')
end

% Workdirecotyr
 pathname = [pwd,'\'];
 % or pathname = 'C:\...\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ IMAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

DP_output = DosePainting;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE CONTOURS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO omvormen naar een functie

% Extract VOI from PET
    PET_VOI = double(PET) .* boundingbox_var;
    save('PET_VOI.mat','PET_VOI')
    save('BoundingBox.mat','boundingbox_var')

% PET max in VOI
    PET_VOI_MAX = max(max(max(PET_VOI)));

% Calculate VOI 50
    % PET 50 threshold
        THRES_PET_50 = PET_VOI_MAX*0.50;

    % PET 50 thresholds    
        PET_VOI_SEG_50 = imbinarize(PET_VOI,THRES_PET_50);

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
% Calculate center of multiple blobs if necessary TODO omvormen naar een
% functie
    
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
        
    % Merge into one variable
    
    centroids = [VOI50_centroid;...
        VOI60_centroid;...
        VOI70_centroid;...
        VOI80_centroid;...
        VOI90_centroid;...
        VOI95_centroid];
    
    
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

% Generate normal vectors
normalvectors = generate_normvect_3arcs(angles,couch_pos);

% Plot normal vectors
if 0
    plot3([0,0],[0,0],[0,1],'k*-',...
        [0,0],[0,0.3640],[0,1],'b^-',...
        [0,0],[0,0.8391],[0,1],'g>-',...
        [0,0],[0,1.7321],[0,1],'r^-',...
        [0,0],[0,5.6713],[0,1],'y>-',...
        [0,0],[0,5.6713],[0,-1],'k^-',...
        [0,0],[0,1.7321],[0,-1],'m>-');
    legend({'0', '20', '40','60','80','100','120'})
    set(gca,'DataAspectRatio',[1 1 1])
end

% Detect "tumour" points for each VOI
original_all = findTumour(PET_VOI_SEG_LARG_50,... 
    PET_VOI_SEG_FILL_60,...
    PET_VOI_SEG_FILL_70,...
    PET_VOI_SEG_FILL_80,...
    PET_VOI_SEG_FILL_90,...
    PET_VOI_SEG_FILL_95);

original_50 = original_all{1};
original_60 = original_all{2};
original_70 = original_all{3};
original_80 = original_all{4};
original_90 = original_all{5};
original_95 = original_all{6};

% Projection
projection_50 = projection(normalvectors,VOI50_centroid,original_50);
projection_60 = projection(normalvectors,VOI60_centroid,original_60);
projection_70 = projection(normalvectors,VOI70_centroid,original_70);
projection_80 = projection(normalvectors,VOI80_centroid,original_80);
projection_90 = projection(normalvectors,VOI90_centroid,original_90);
projection_95 = projection(normalvectors,VOI95_centroid,original_95);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT PROJECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    HAxes = plotProjection(normalvectors(4,:),VOI50_centroid,original_50,projection_50{4});
end 
%TODO use view option for better results

if 0    % Addition of bounding box points; adjust 
    %11 G60 C-45
    plot3(229.1591, 197.5494, 111.2419,'hc','Parent',HAxes)
    plot3(229.1591, 212.4065, 85.5086,'hc','Parent',HAxes)
    plot3(244.9294, 197.5494, 83.9269,'hc','Parent',HAxes)
    plot3(244.9294, 212.4065, 58.1936,'hc','Parent',HAxes)  

    % 1 G0 C0
    plot3(228, 196, 83.7518,'hc','Parent',HAxes) % cyan
    plot3(228, 215, 83.7518,'hm','Parent',HAxes) % magenta
    plot3(248, 196, 83.7518,'hg','Parent',HAxes) % green
    plot3(248, 215, 83.7518,'hy','Parent',HAxes) % yellow 
    % 4 G60 C0
    plot3(228, 200.8522, 91.4920,'hg','Parent',HAxes) % cyan
    plot3(228, 210.0973, 75.4789,'hg','Parent',HAxes) % magenta
    plot3(248, 200.8522, 91.4920,'hg','Parent',HAxes) % green
    plot3(248, 210.0973, 75.4789,'hg','Parent',HAxes) % yellow
    %8 G0 C-45
    plot3(228, 196, 83.7518,'hc','Parent',HAxes) % cyan
    plot3(228, 215, 83.7518,'hm','Parent',HAxes) % magenta
    plot3(248, 196, 83.7518,'hg','Parent',HAxes) % green
    plot3(248, 215, 83.7518,'hy','Parent',HAxes) % yellow   
    %11 G60 C-45
    plot3(250.1762, 197.5337, 74.8030,'hc','Parent',HAxes) % cyan
    plot3(235.3190, 212.3908, 74.8030,'hm','Parent',HAxes) % magenta
    plot3(239.5645, 197.5337, 93.1830,'hg','Parent',HAxes) % green
    plot3(224.7073, 212.3908, 93.1830,'hb','Parent',HAxes) % yellow      
    % 15 G0 C-90
    plot3(228, 196, 83.7518,'hc','Parent',HAxes) % cyan
    plot3(228, 215, 83.7518,'hm','Parent',HAxes) % magenta
    plot3(248, 196, 83.7518,'hg','Parent',HAxes) % green
    plot3(248, 215, 83.7518,'hy','Parent',HAxes) % yellow
    % 18 C60 C-90
    plot3(231.8854, 196, 93.0886,'hc','Parent',HAxes) % cyan
    plot3(231.8854, 215, 93.0886,'hm','Parent',HAxes) % magenta
    plot3(242.1127, 196, 75.3745,'hg','Parent',HAxes) % green
    plot3(242.1127, 215, 75.3745,'hy','Parent',HAxes) % yellow    
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BOUNDING BOX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define variables
SROW = [PET_info.raw.srow_x; ...
    PET_info.raw.srow_y; ...
    PET_info.raw.srow_z; ...
    0 0 0 1]; % Transformation matrix

% Generate boundaries
% Input (volgorde) projections belangrijk!!
output = jawPos_3arcs(angles,couch_pos,centroids,SROW,...
    projection_50,...
    projection_60,...
    projection_70,...
    projection_80,...
    projection_90,...
    projection_95);

save('SARRPinput.mat','output')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INFORMATION TABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gantry_angle = zeros(length(angles)*6,1);       % 6 different SUV levels
Couch_position = zeros(length(angles)*6,1);
VOI = cell(length(angles)*6,1);
Yawdist_X = zeros(length(angles)*6,1);
Yawdist_Y = zeros(length(angles)*6,1);
Yaw_Mid = zeros(length(angles)*6,3);
counter_table = 1;

fn_VOI = fieldnames(output.beam);
for k = 1:length(angles)
   for m = 4:length(fn_VOI)
                Gantry_angle(counter_table)     = output.beam(k).angle;
                Couch_position(counter_table)   = output.beam(k).couch_pos;
                VOI{counter_table}              = fn_VOI{m};
                Yawdist_X(counter_table)        = output.beam(k).(fn_VOI{m}).('YawDist_X');
                Yawdist_Y(counter_table)        = output.beam(k).(fn_VOI{m}).('YawDist_Y');
                Yaw_Mid(counter_table,:)        = output.beam(k).(fn_VOI{m}).('YawMid_SARRP');
                counter_table       = counter_table + 1;
   end
end

SARRPinput = table(Gantry_angle,Couch_position,VOI,Yawdist_X,Yawdist_Y,Yaw_Mid)
