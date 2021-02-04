function [projected_points] = orthoProject_v2(normvect,point,points_to_project)
% input variables:
%   - Variables plane:
%           - Normal vector: (a,b,c)
%           - Point: (x1,y1,z1)
%   - point to project: (x0,y0,z0)
%
% output variables:
%   - projected point: [xp, yp, zp]
%TODO matrix multiplication for efficiency?

% Defining variables
a = normvect(1);
b = normvect(2);
c = normvect(3);
x_1 = point(1);
y_1 = point(2);
z_1 = point(3);
n = size(points_to_project,1);  % number of points to project
projected_points = zeros(n,3);

for i = 1:n
    x_0 = points_to_project(i,1);
    y_0 = points_to_project(i,2);
    z_0 = points_to_project(i,3);
    
    t = (-a*x_0 + a*x_1 - b*y_0 + b*y_1 - c*z_0 + c*z_1)/(a^2 + b^2 + c^2);

    projected_points(i,1) = x_0 + a * t;
    projected_points(i,2) = y_0 + b * t;
    projected_points(i,3) = z_0 + c * t;
end
end