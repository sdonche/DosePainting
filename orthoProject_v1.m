function [x_p,y_p,z_p] = orthoProject_v1(a,b,c,x_1,y_1,z_1,x_0,y_0,z_0)
% input variables:
%   - Variables plane:
%           - Normal vector: (a,b,c)
%           - Point: (x1,y1,z1)
%   - point to project: (x0,y0,z0)
%
% output variables:
%   - projected point: [xp, yp, zp]

t = (-a*x_0 + a*x_1 - b*y_0 + b*y_1 - c*z_0 + c*z_1)/(a^2 + b^2 + c^2);

x_p = x_0 + a * t;
y_p = y_0 + b * t;
z_p = z_0 + c * t;
end