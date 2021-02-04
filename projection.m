function [projected_points] = projection(normvectors,centroid,points_to_project)
%TODO UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Defining parameters
num_vect = size(normvectors,1);
projected_points = cell(num_vect,1);

for i = 1:num_vect

    projected_points{i} = orthoProject_v2(normvectors(i,:),centroid,points_to_project);

end
end