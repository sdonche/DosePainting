function NewImageMask = ExtractSegments(ImageMask, numberOfSegments)
% ExtractSegment extracts the desired segment from the mask image.
% The absolute value of numberOfSegments gives the amount of segments that
% will be extracted.
% If numberOfSegments is positive (> 0) then the largest segment(s) are
% extracted.
% If numberOfSegments is negatieve (< 0) then the smallest segment(s) are
% extracted.
% A mask with the desired number of segments is returned.
%try
  % Extract segment properties  
  [labeledImage, numberOfBlobs] = bwlabeln(ImageMask);
  blobMeasurements = regionprops3(labeledImage, 'volume');
  % Get all the areas
  allVolumes = [blobMeasurements.Volume];
  
  if numberOfSegments > 0
    % For positive numbers, sort in order of largest to smallest
    [sortedAreas, sortIndexes] = sort(allVolumes, 'descend');
  elseif numberOfSegments < 0
    % For negative numbers, sort in order of smallest to largest
    [sortedAreas, sortIndexes] = sort(allVolumes, 'ascend');
    % Need to negate numberToExtract for later
    numberToExtract = -numberToExtract;
  else
    % numberToExtract = 0. Return no segments
    NewImageMask = false(size(ImageMask));
    return;
  end
  % Extract the "numberToExtract" largest blob(a)s using ismember()
  BiggestSegment = ismember(labeledImage, sortIndexes(1:numberOfSegments));
  % Convert from integer labeled image into binary (logical) image.
  NewImageMask = BiggestSegment > 0;
%catch ME
%   errorMessage = sprintf('Error in function ExtractNLargestBlobs().\n\nError Message:\n%s', ME.message);
%   fprintf(1, '%s\n', errorMessage);
%   uiwait(warndlg(errorMessage));
end
