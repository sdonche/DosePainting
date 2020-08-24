%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Dose Painting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Information
% Determination SARRP isocenters
% Determined isocenters:
%           Isocenter 1: 
%                   MRI-based
%                   A box of 10 mm x 10 mm x 10 mm is dragged around the
%                   T1ce region to fit as good as possible. Isocenter 1
%                   corresponds to the middlepoint of this box.
%           Isocenter 2:
%                   PET-based
%                   A volume is determined to correspond to 5% highest SUV
%                   uptake values. Isocenter 2 corresponds with the middle
%                   point of this volume.
%
% Conform PhD Sam Donche

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLEAN SLATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
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

disp('Reading Toolboxes...')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Extra\spm12')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Extra\NIfTI_20140122')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Extra\Medical Image Reader and Viewer')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Extra\nifti_utils-master\nifti_utils') % "updated" version of NIfTI_20140122
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\imshow3D.m')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\imshow3Dfull.m')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Preclinical\Scripts\App_BoundingBox')
addpath('C:\Users\Hoofdgebruiker\OneDrive - UGent\Doctoraat\MATLAB\Preclinical\Scripts')

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
%% SET BOUNDING BOX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The MRI is used as a bounding box around the tumour

application = DosePainting;

%(MRI); % Don't close application!
%application2 = BoundingBox(MRI)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MRI ISOCENTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MRI ceT1 region: 0.05/0.05% highest percentile?

mri_roi = double(MRI) .* boundingbox_var;
mri_roi_nonzero = mri_roi > 0;

threshold = prctile(mri_roi(mri_roi_nonzero),95);

% mri_roi_max = max(max(max(mri_roi)));
% threshold = mri_roi_max * 0.95;

%mri_roi_segment = mri_roi > threshold;
mri_roi_segment = imbinarize(mri_roi,threshold);

% figure()
% sliceViewer(mri_roi)
% figure()
% sliceViewer(mri_roi_segment)

% Fill holes
mri_roi_segment_filled = imfill(mri_roi_segment, 'holes');

% figure()
% sliceViewer(mri_roi_segment)
% figure()
% sliceViewer(mri_roi_segment_filled)

% Extract properties
[labeledImage, numberOfBlobs] = bwlabeln(mri_roi_segment_filled);
blobMeasurements_MRI = regionprops3(labeledImage, 'Volume');

% Extract largest segment
mri_roi_segment_largest = ExtractSegments(mri_roi_segment_filled,1);

figure()
sliceViewer(mri_roi_segment_filled)
figure()
sliceViewer(mri_roi_segment_largest)

% Extract center of largest segment

[labeledImage, numberOfBlobs] = bwlabeln(mri_roi_segment_largest); % TODO check if numberOfBlobs is 1
blobMeasurements_MRI = regionprops3(labeledImage, mri_roi, 'Volume', 'Centroid', 'WeightedCentroid', ... 
    'MaxIntensity','MeanIntensity','MinIntensity','Orientation','BoundingBox','EquivDiameter'); 

blobMeasurements_MRI.Centroid % TODO see what kind of center this is

% Show point on MRI an CT
% Center volume
xslice = blobMeasurements_MRI.Centroid(1);
yslice = blobMeasurements_MRI.Centroid(2);
zslice = blobMeasurements_MRI.Centroid(3);
% MRI
% 
figure('units','normalized','outerposition',[0 0 1 1])
subplot(3,3,1)
imshow(squeeze(MRI(:,round(xslice),:)),[0 25000]);
view(90,90)
axis on
hold on
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(2),'o')
axis off
hold off
% 
subplot(3,3,2)
imshow(squeeze(MRI(round(yslice),:,:)),[0 25000]);
view(90,90)
axis on
hold on
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(1),'o')
axis off
hold off
% XY-plane
subplot(3,3,3)
imshow(squeeze(MRI(:,:,round(zslice))),[0 25000]);
axis on
hold on
plot(blobMeasurements_MRI.Centroid(1),blobMeasurements_MRI.Centroid(2),'o')
axis off
hold off
%CT
CT_max = max(max(max(CT)));
%
subplot(3,3,4)
imshow(squeeze(CT(:,round(xslice),:)),[0.35*CT_max 0.50*CT_max]);
view(90,90)
axis on
hold on
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(2),'o')
axis off
hold off
% 
subplot(3,3,5)
imshow(squeeze(CT(round(yslice),:,:)),[0.35*CT_max 0.50*CT_max]);
view(90,90)
axis on
hold on
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(1),'o')
axis off
hold off
% XY-plane
subplot(3,3,6)
imshow(squeeze(CT(:,:,round(zslice))),[0.35*CT_max 0.50*CT_max]);
axis on
hold on
plot(blobMeasurements_MRI.Centroid(1),blobMeasurements_MRI.Centroid(2),'o')
axis off
hold off
%PET
%
subplot(3,3,7)
imshow(squeeze(PET(:,round(xslice),:)));
view(90,90)
axis on
hold on
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(2),'o')
axis off
hold off
% 
subplot(3,3,8)
imshow(squeeze(PET(round(yslice),:,:)));
view(90,90)
axis on
hold on
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(1),'o')
axis off
hold off
% XY-plane
subplot(3,3,9)
imshow(squeeze(PET(:,:,round(zslice))));
axis on
hold on
plot(blobMeasurements_MRI.Centroid(1),blobMeasurements_MRI.Centroid(2),'o')
axis off
hold off




figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
imshowpair(squeeze(CT(:,:,round(zslice))),squeeze(MRI(:,:,round(zslice))))
axis on
hold on
plot(blobMeasurements_MRI.Centroid(1),blobMeasurements_MRI.Centroid(2),'o')
axis off
hold off
subplot(1,2,2)
imshowpair(squeeze(MRI(:,round(xslice),:)),squeeze(CT(:,round(xslice),:)))
view(90,90)
axis on
hold on
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(2),'o')
axis off
hold off

% TODO: plaats 10x10x10 box rond punt en toon

% Convert pixel center to world center

SROW = [CT_info.raw.srow_x; CT_info.raw.srow_y; CT_info.raw.srow_z; 0 0 0 1];
P_COOR = [blobMeasurements_MRI.Centroid(1); blobMeasurements_MRI.Centroid(2); blobMeasurements_MRI.Centroid(3); 1];

R_COOR = SROW * P_COOR;

% Convert world coordinates to Muriplan coordinates

% MURI_COOR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PET ISOCENTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO:  'ScrollWheelDuringDraw' misschien handig voor tekenen ROI in app?

% Info from MRI
% Center = blobMeasurements_MRI.Centroid;
% Diameter = blobMeasurements_MRI.EquivDiameter;
% Diameter_PET = Diameter*1.5;
% Radius_PET = Diameter_PET/2;

% Volume MRI ROI expanderen
se = strel('sphere',10); %TODO eventueel aanpassen
dilated_mask = imdilate(mri_roi_segment_largest, se);

% PET ROI
pet_roi = double(PET) .* double(dilated_mask);

pet_max = max(max(max(PET)));
pet_roi_max = max(max(max(pet_roi)));


% figure()
% sliceViewer(PET,'DisplayRange', [0 pet_max])
% figure()
% sliceViewer(pet_roi,'DisplayRange', [0 pet_max])

% pet_roi_nonzero = pet_roi > 0;
% threshold_pet2 = prctile(pet_roi(pet_roi_nonzero)) incorrect
threshold_pet = pet_roi_max*0.95;

%pet_roi_segment = pet_roi > threshold_pet;
pet_roi_segment = imbinarize(pet_roi,threshold_pet);

% figure()
% sliceViewer(pet_roi)
% figure()
% sliceViewer(pet_roi_segment,'DisplayRange',[0 1])


% Fill holes
pet_roi_segment_filled = imfill(pet_roi_segment, 'holes');

% figure()
% sliceViewer(pet_roi_segment)
% figure()
% sliceViewer(pet_roi_segment_filled)

% Extract properties
[labeledImage, numberOfBlobs] = bwlabeln(pet_roi_segment_filled);
blobMeasurements_PET = regionprops3(labeledImage, 'Volume');

% Extract largest segment
pet_roi_segment_largest = ExtractSegments(pet_roi_segment_filled,1);

% figure()
% sliceViewer(mri_roi_segment_filled)
% figure()
% sliceViewer(mri_roi_segment_largest)

% Extract center of largest segment

[labeledImage, numberOfBlobs] = bwlabeln(pet_roi_segment_largest); % TODO check if numberOfBlobs is 1
blobMeasurements_PET = regionprops3(labeledImage, pet_roi, 'Volume', 'Centroid', 'WeightedCentroid', ... 
    'MaxIntensity','MeanIntensity','MinIntensity','Orientation','BoundingBox','EquivDiameter'); 

blobMeasurements_PET.Centroid

% Show point on MRI an CT
% Center volume
xslice = blobMeasurements_PET.Centroid(1);
yslice = blobMeasurements_PET.Centroid(2);
zslice = blobMeasurements_PET.Centroid(3);
% MRI
% 
figure()
subplot(3,3,1)
imshow(squeeze(MRI(:,round(xslice),:)),[0 25000]);
view(90,90)
axis on
hold on
plot(blobMeasurements_PET.Centroid(3),blobMeasurements_PET.Centroid(2),'or')
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(2),'ob')
axis off
hold off
% 
subplot(3,3,2)
imshow(squeeze(MRI(round(yslice),:,:)),[0 25000]);
view(90,90)
axis on
hold on
plot(blobMeasurements_PET.Centroid(3),blobMeasurements_PET.Centroid(1),'or')
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(1),'ob')
axis off
hold off
% XY-plane
subplot(3,3,3)
imshow(squeeze(MRI(:,:,round(zslice))),[0 25000]);
axis on
hold on
plot(blobMeasurements_PET.Centroid(1),blobMeasurements_PET.Centroid(2),'or')
plot(blobMeasurements_MRI.Centroid(1),blobMeasurements_MRI.Centroid(2),'ob')
axis off
hold off
%CT
CT_max = max(max(max(CT)));
%
subplot(3,3,4)
imshow(squeeze(CT(:,round(xslice),:)),[0.35*CT_max 0.50*CT_max]);
view(90,90)
axis on
hold on
plot(blobMeasurements_PET.Centroid(3),blobMeasurements_PET.Centroid(2),'or')
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(2),'ob')
axis off
hold off
% 
subplot(3,3,5)
imshow(squeeze(CT(round(yslice),:,:)),[0.35*CT_max 0.50*CT_max]);
view(90,90)
axis on
hold on
plot(blobMeasurements_PET.Centroid(3),blobMeasurements_PET.Centroid(1),'or')
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(1),'ob')
axis off
hold off
% XY-plane
subplot(3,3,6)
imshow(squeeze(CT(:,:,round(zslice))),[0.35*CT_max 0.50*CT_max]);
axis on
hold on
plot(blobMeasurements_PET.Centroid(1),blobMeasurements_PET.Centroid(2),'or')
plot(blobMeasurements_MRI.Centroid(1),blobMeasurements_MRI.Centroid(2),'ob')
axis off
hold off
%PET
%
subplot(3,3,7)
imshow(squeeze(PET(:,round(xslice),:)));
view(90,90)
axis on
hold on
plot(blobMeasurements_PET.Centroid(3),blobMeasurements_PET.Centroid(2),'or')
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(2),'ob')
axis off
hold off
% 
subplot(3,3,8)
imshow(squeeze(PET(round(yslice),:,:)));
view(90,90)
axis on
hold on
plot(blobMeasurements_PET.Centroid(3),blobMeasurements_PET.Centroid(1),'or')
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(1),'ob')
axis off
hold off
% XY-plane
subplot(3,3,9)
imshow(squeeze(PET(:,:,round(zslice))));
axis on
hold on
plot(blobMeasurements_PET.Centroid(1),blobMeasurements_PET.Centroid(2),'or')
plot(blobMeasurements_MRI.Centroid(1),blobMeasurements_MRI.Centroid(2),'ob')
axis off
hold off




figure()
subplot(1,2,1)
imshowpair(squeeze(CT(:,:,round(zslice))),squeeze(MRI(:,:,round(zslice))))
axis on
hold on
plot(blobMeasurements_PET.Centroid(1),blobMeasurements_PET.Centroid(2),'or')
plot(blobMeasurements_MRI.Centroid(1),blobMeasurements_MRI.Centroid(2),'ob')
axis off
hold off
subplot(1,2,2)
imshowpair(squeeze(MRI(:,round(xslice),:)),squeeze(CT(:,round(xslice),:)))
view(90,90)
axis on
hold on
plot(blobMeasurements_PET.Centroid(3),blobMeasurements_PET.Centroid(2),'or')
plot(blobMeasurements_MRI.Centroid(3),blobMeasurements_MRI.Centroid(2),'ob')
axis off
hold off

% TODO: plaats 10x10x10 box rond punt en toon

% Convert pixel center to world center

SROW = [CT_info.raw.srow_x; CT_info.raw.srow_y; CT_info.raw.srow_z; 0 0 0 1];
P_COOR = [blobMeasurements_PET.Centroid(1); blobMeasurements_PET.Centroid(2); blobMeasurements_PET.Centroid(3); 1];

R_COOR = SROW * P_COOR;

% Convert world coordinates to Muriplan coordinates

% MURI_COOR




%PET
h = sliceViewer(PET)
hold on
%axes on
i = drawcuboid('Position', [Center(1)-Radius_PET,Center(2)-Radius_PET,Center(3)-Radius_PET,Diameter_PET,Diameter_PET,Diameter_PET], ...
    'InteractionsAllowed', 'none', 'Label','Bounding Box PET','LabelVisible','hover');

g = imshow(squeeze(PET(:,233,:)))
hold on
%axes on
view(90,90)
j = drawrectangle('Position', [Center(1)-Radius_PET,Center(3)-Radius_PET,Diameter_PET,Diameter_PET], ...
    'InteractionsAllowed', 'none', 'Label','Bounding Box PET','LabelVisible','hover');




% TODO: if extractlagergest segment fails -> extract ## largest segments
% and delete the appropriate segments

% TODO : 
% 'WeightedCentroid'	
% Center of the region based on location and intensity value, 
% returned as a p-by-3 vector of coordinates. 
% The first element of WeightedCentroid is the horizontal coordinate (or x-coordinate) 
% of the weighted centroid. The second element is the vertical coordinate (or y-coordinate). 
% The third element is the planar coordinate (or z-coordinate).
% 'BoundingBox'	
% Smallest cuboid containing the region, returned as a 1-by-6 vector of the form 
% [ulf_x ulf_y ulf_z width_x width_y width_z]. ulf_x, ulf_y, and ulf_z 
% specify the upper-left front corner of the cuboid. width_x, width_y, 
% and width_z specify the width of the cuboid along each dimension.
% 'Centroid'	
% Center of mass of the region, returned as a 1-by-3 vector of the form 
% [centroid_x centroid_y and centroid_z]. The first element, centroid_x, 
% is the horizontal coordinate (or x-coordinate) of the center of mass. 
%The second element, centroid_y, is the vertical coordinate (or y-coordinate). 
% The third element, centroid_z, is the planar coordinate (or z-coordinate).



% Probeersels hieronder


[counts,edges] = histcounts(mri_roi(mri_roi_nonzero));

[threshold, effectiveness] = otsuthresh(counts);


edges = 0:15:3000;
histogram(mri_roi(mri_roi_nonzero), edges);

h = histogram(mri_roi(mri_roi_nonzero));

h.BinEdges


figure()
bar(counts)



% Place vertical bar on histogram to show threshold
yl = ylim();
line([threshold, threshold], [yl(1), yl(2)], ...
  'Color', 'r', 'LineWidth', 2);
grid on;

    % Read an image of coins and compute a 16-bin histogram.
    I = imread('coins.png');
    counts = imhist(I, 16);
 
    % Compute a global threshold using the histogram counts.
    T = otsuthresh(counts);
 
    % Binarize image using computed threshold.
    BW = imbinarize(I,T);
 
    figure, imshow(BW)






% Display it
figure; % new figure
subplot(2, 2, 1);
bar(diffHistogram, 'BarWidth', 1);
title('Difference of Histograms', 'FontSize', fontSize);
% Now let's do an Otsu thresholding on it.
thresholdLevel = 255 * graythresh(diffHistogram)  % Find Otsu threshold level
% Place vertical bar on histogram to show threshold
yl = ylim();
line([thresholdLevel, thresholdLevel], [yl(1), yl(2)], ...
  'Color', 'r', 'LineWidth', 2);
grid on;
% Find pixels above the threshold
aboveThreshold = grayImage > thresholdLevel;
subplot(2, 2, 2);
imshow(aboveThreshold, []);
title('Above threshold', 'FontSize', fontSize);



figure()
g = sliceViewer(MRI,'DisplayRange',[0 1000])
h = drawcuboid(gca)


PETbox = PET;

for i = 1:length(Boxx)
    PETbox(Boxx(i),Boxy(i),Boxz(i)) = 50000; 
end

CTbox = CT;
for i = 1:length(Boxx)
    CTbox(Boxx(i),Boxy(i),Boxz(i)) = 50000; 
end

figure()
sliceViewer(PETbox)

figure()
sliceViewer(CTbox)

figure()
sliceViewer(MRI)




% hFig = figure;  % Bring up a new figure.
% for k = 1 : size(MRI, 3)
%   subplot(1, 3, 1);
%   cla; % clear current axis
%   imshow(MRI(:,:,i), []);
%   drawnow;
%   if k == 1
%     hFH = drawfreehand(); % imfreehand
%     binaryImage = hFH.createMask;
%     subplot(1, 3, 2);
%     cla;
%     imshow(binaryImage);
%     % Enlarge figure to full screen.
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     drawnow;
%     % xy = hFH.getPosition();  % Unused
%   end
%   blackMaskedImage = MRI(:,:,i) .* binaryImage(:,:,i);
%   subplot(1, 3, 3);
%   cla;
%   imshow(blackMaskedImage);
%   drawnow;
% end

% You can specify a region of interest in many ways.

% You can interactively draw it using tools like 
% roipoly(): select polygon region of interest within image
% roipolyold(): old version of above
% rbbox(): 
% imrect(): use drawrectangle instead; draggable rectangle
% imfreehand()
% You could specify a box by hard coding in the rows and columns. 
% You could do segmentation, for example by intensity thresholding or color segmentation.
% However you do it, you have a binary image that is the mask. You can get a frame of your video 
% and then mask it with the ROI mask like this:

rgbImage = MRI(:,:,1);
mask = zeros(size(MRI(:,:,1)));

% Mask the image.
maskedRgbImage = bsxfun(@times, rgbImage, cast(mask,class(rgbImage)));

% Display image
subplot(2,2,1);
imshow(rgbImage);
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Have user freehand draw a mask in the image.
uiwait(msgbox('Draw a region with left button down.  Double click inside to finish it'));
h = imfreehand();
vertices = wait(h);
% Create a binary mask
mask = createMask(h);
subplot(2,2,2);
imshow(mask, []);
% Mask the image.
maskedRgbImage = bsxfun(@times, rgbImage, cast(mask,class(rgbImage)));
subplot(2,2,3);
imshow(maskedRgbImage);


f = figure();
imshow(MRI(:,:,1))
a = axes('Parent',f);
cla;
k = waitforbuttonpress;
point1 = a.CurrentPoint;    % button down detected
finalRect = rbbox;                 % return figure units
point2 = a.CurrentPoint;    % button up detected
point1 = point1(1,1:2);            % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);           % calculate locations
offset = abs(point1-point2);       % and dimensions
x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
hold on
axis manual
plot(x,y)                          % redraw in dataspace units




figure, imshow3Dfull(MRI,[0 15000]);
h = imrect(gca, [10 10 100 100]);
addNewPositionCallback(h,@(p) title(mat2str(p,3)));
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(h,fcn);




[slice_1_roi tmp_xi tmp_yi]=roipoly(squeeze(MRI(:,48,:)));  
imshow(squeeze(MRI(:,48,:))); %To get the figure handle.
h = impoly(gca,[tmp_xi tmp_yi]); %Alter your mask here. DONT close the figure window yet.
BW = createMask(h); %BW contains the mask which you just altered.
pos = getPosition(h); %Stores the points in your new mask.
%Close the figure now, and repeat this for the next slice.
imshow(D(:,:,3));
h = impoly(gca,pos); %Now pos contains the points from your mask.



image_test = squeeze(MRI(:,234,:));

figure()
h = imshow(image_test);
g = drawrectangle('Label','Bounding Box','Position',[10 10 54.5455 54.5455],'InteractionsAllowed','translate');
mask = createMask(g);

newimage = image_test .* int16(mask);

figure()
imshow(newimage);



grayImage = squeeze(MRI(:,234,:));

% Get the dimensions of the image.
% numberOfColorBands should be = 1.
[rows, columns, numberOfColorBands] = size(grayImage);
% if numberOfColorBands > 1
%   grayImage = rgb2gray(grayImage);
% end
% Display the original gray scale image.
subplot(2, 3, 1);
imshow(grayImage, []);
title('Original Grayscale Image', 'FontSize', 16);
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% Give a name to the title bar.
set(gcf,'name','Demo by ImageAnalyst','numbertitle','off')
% Let's compute and display the histogram.

subplot(2, 3, 4);
h = histogram(grayImage,1000);
grid on;
title('Histogram of original image', 'FontSize', 16);

% [pixelCount, grayLevels] = imhist(grayImage, 1000);
% subplot(2, 3, 4);
% bar(pixelCount);
% grid on;
% title('Histogram of original image', 'FontSize', 16);
% xlim([0 median(grayLevels)]); % Scale x axis manually.

% Extract left half and display histogram
middleColumn = floor(columns/2);
leftHalfImage = grayImage(:, 1:middleColumn);
subplot(2, 3, 2);
imshow(leftHalfImage, []);
title('left Half Image', 'FontSize', 16);
% Let's compute and display the histogram.
[pixelCountL, grayLevelsL] = imhist(leftHalfImage, 256);
subplot(2, 3, 5);
bar(pixelCountL);
grid on;
title('Histogram of left half image', 'FontSize', fontSize);
xlim([0 grayLevels(end)]); % Scale x axis manually.
% Extract right half and display histogram
middleColumn = floor(columns/2);
rightHalfImage = grayImage(:, middleColumn+1:end);
subplot(2, 3, 3);
imshow(rightHalfImage, []);
title('Right Half Image', 'FontSize', fontSize);
% Let's compute and display the histogram.
[pixelCountR, grayLevelsR] = imhist(rightHalfImage, 256);
subplot(2, 3, 6);
bar(pixelCountR);
grid on;
title('Histogram of right half image', 'FontSize', fontSize);
xlim([0 grayLevels(end)]); % Scale x axis manually.
% Now let's subtract the histograms (for some reason);
diffHistogram = int16(pixelCountL - pixelCountR);
% Display it
figure; % new figure
subplot(2, 2, 1);
bar(diffHistogram, 'BarWidth', 1);
title('Difference of Histograms', 'FontSize', fontSize);
% Now let's do an Otsu thresholding on it.
thresholdLevel = 255 * graythresh(diffHistogram)  % Find Otsu threshold level
% Place vertical bar on histogram to show threshold
yl = ylim();
line([thresholdLevel, thresholdLevel], [yl(1), yl(2)], ...
  'Color', 'r', 'LineWidth', 2);
grid on;
% Find pixels above the threshold
aboveThreshold = grayImage > thresholdLevel;
subplot(2, 2, 2);
imshow(aboveThreshold, []);
title('Above threshold', 'FontSize', fontSize);


I = imread('coins.png');
subplot(2,1,1)
imshow(I)
level = graythresh(I);
BW = imbinarize(I,level);
subplot(2,1,2)
imshow(BW)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract biggest blob
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function ExtractBiggestBlob()
% clc;    % Clear the command window.
% close all;  % Close all figures (except those of imtool.)
% imtool close all;  % Close all imtool figures if you have the Image Processing Toolbox.
% clear;  % Erase all existing variables. Or clearvars if you want.
% workspace;  % Make sure the workspace panel is showing.
% format long g;
% format compact;
% fontSize = 20;
% % Read in a standard MATLAB gray scale demo image.
% folder = fullfile(matlabroot, '\toolbox\images\imdemos');
% baseFileName = 'coins.png';
% % Get the full filename, with path prepended.
% fullFileName = fullfile(folder, baseFileName);
% % Check if file exists.
% if ~exist(fullFileName, 'file')
%   % File doesn't exist -- didn't find it there.  Check the search path for it.
%   fullFileName = baseFileName; % No path this time.
%   if ~exist(fullFileName, 'file')
%     % Still didn't find it.  Alert user.
%     errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
%     uiwait(warndlg(errorMessage));
%     return;
%   end
% end
% grayImage = imread(fullFileName);
% % Get the dimensions of the image.
% % numberOfColorBands should be = 1.
% [rows, columns, numberOfColorBands] = size(grayImage);
% % Display the original gray scale image.
% subplot(2, 2, 1);
% imshow(grayImage, []);
% title('Original Grayscale Image', 'FontSize', fontSize);
% % Enlarge figure to full screen.
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% % Give a name to the title bar.
% set(gcf,'name','Demo by ImageAnalyst','numbertitle','off')
% % Let's compute and display the histogram.
% [pixelCount, grayLevels] = imhist(grayImage);
% subplot(2, 2, 2);
% bar(pixelCount);
% grid on;
% title('Histogram of original image', 'FontSize', fontSize);
% xlim([0 grayLevels(end)]); % Scale x axis manually.
% % Threshold the image to binarize it.
% binaryImage = grayImage > 100;
% % Fill holes
% binaryImage = imfill(binaryImage, 'holes');
% % Display the image.
% subplot(2, 2, 3);
% imshow(binaryImage, []);
% title('Binary Image', 'FontSize', fontSize);
% % Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
% [labeledImage, numberOfBlobs] = bwlabel(binaryImage); % 3d variant: bwlabeln
% blobMeasurements = regionprops(labeledImage, 'area', 'Centroid'); % 3d variant: regionprops3
% % Get all the areas
% allAreas = [blobMeasurements.Area] % No semicolon so it will print to the command window.
% menuOptions{1} = '0'; % Add option to extract no blobs.
% % Display areas on image
% for k = 1 : numberOfBlobs           % Loop through all blobs.
%   thisCentroid = [blobMeasurements(k).Centroid(1), blobMeasurements(k).Centroid(2)];
%   message = sprintf('Area = %d', allAreas(k));
%   text(thisCentroid(1), thisCentroid(2), message, 'Color', 'r');
%   menuOptions{k+1} = sprintf('%d', k);
% end
% % Ask user how many blobs to extract.
% numberToExtract = menu('How many do you want to extract', menuOptions) - 1;
% % Ask user if they want the smallest or largest blobs.
% promptMessage = sprintf('Do you want the %d largest, or %d smallest, blobs?',...
%   numberToExtract, numberToExtract);
% titleBarCaption = 'Largest or Smallest?';
% sizeOption = questdlg(promptMessage, titleBarCaption, 'Largest', 'Smallest', 'Cancel', 'Largest');
% if strcmpi(sizeOption, 'Cancel')
%   return;
% elseif strcmpi(sizeOption, 'Smallest')
%   % If they want the smallest, make the number negative.
%   numberToExtract = -numberToExtract;
% end
% %---------------------------------------------------------------------------
% % Extract the largest area using our custom function ExtractNLargestBlobs().
% % This is the meat of the demo!
% biggestBlob = ExtractNLargestBlobs(binaryImage, numberToExtract);
% %---------------------------------------------------------------------------
% % Display the image.
% subplot(2, 2, 4);
% imshow(biggestBlob, []);
% % Make the number positive again.  We don't need it negative for smallest extraction anymore.
% numberToExtract = abs(numberToExtract);
% if numberToExtract == 1
%   caption = sprintf('Extracted %s Blob', sizeOption);
% elseif numberToExtract > 1
%   caption = sprintf('Extracted %d %s Blobs', numberToExtract, sizeOption);
% else % It's zero
%   caption = sprintf('Extracted 0 Blobs.');
% end
% title(caption, 'FontSize', fontSize);
% msgbox('Done with demo!');
% 
% 
% % Function to return the specified number of largest or smallest blobs in a binary image.
% % If numberToExtract > 0 it returns the numberToExtract largest blobs.
% % If numberToExtract < 0 it returns the numberToExtract smallest blobs.
% % Example: return a binary image with only the largest blob:
% %   binaryImage = ExtractNLargestBlobs(binaryImage, 1)
% % Example: return a binary image with the 3 smallest blobs:
% %   binaryImage = ExtractNLargestBlobs(binaryImage, -3)
% 
% 
% function binaryImage = ExtractNLargestBlobs(binaryImage, numberToExtract)
% try
%   % Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
%   [labeledImage, numberOfBlobs] = bwlabel(binaryImage);
%   blobMeasurements = regionprops(labeledImage, 'area');
%   % Get all the areas
%   allAreas = [blobMeasurements.Area];
%   if numberToExtract > 0
%     % For positive numbers, sort in order of largest to smallest.
%     % Sort them.
%     [sortedAreas, sortIndexes] = sort(allAreas, 'descend');
%   elseif numberToExtract < 0
%     % For negative numbers, sort in order of smallest to largest.
%     % Sort them.
%     [sortedAreas, sortIndexes] = sort(allAreas, 'ascend');
%     % Need to negate numberToExtract so we can use it in sortIndexes later.
%     numberToExtract = -numberToExtract;
%   else
%     % numberToExtract = 0.  Shouldn't happen.  Return no blobs.
%     binaryImage = false(size(binaryImage));
%     return;
%   end
%   % Extract the "numberToExtract" largest blob(a)s using ismember().
%   biggestBlob = ismember(labeledImage, sortIndexes(1:numberToExtract));
%   % Convert from integer labeled image into binary (logical) image.
%   binaryImage = biggestBlob > 0;
% catch ME
%   errorMessage = sprintf('Error in function ExtractNLargestBlobs().\n\nError Message:\n%s', ME.message);
%   fprintf(1, '%s\n', errorMessage);
%   uiwait(warndlg(errorMessage));
% end

