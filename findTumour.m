function [output] = findTumour(varargin)
%TODO UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

n = nargin;
output = cell(1,n);

for i = 1:n
   
    index = find(varargin{i});
    [original_y, original_x, original_z] = ind2sub((size(varargin{i})),index); %ind2sub output = [ROW, COLUMN]!!!
    original = [original_x, original_y, original_z];
    output{i} = original;
    
end
end

