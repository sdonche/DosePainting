function distance = point2line(point,vertex1,vertex2)
%POINT2LINE Calculates distance between point and line
%   "Point" should be nx3
%   "vertex1" and "vertex2" are verteces on the line (each 1x3)
%   "distance" is an nx1 vetor with the orthogonal distances

vertex1 = repmat(vertex1,size(point,1),1);
vertex2 = repmat(vertex2,size(point,1),1);
A = vertex1 - vertex2;
B = point - vertex2;
distance = sqrt(sum(cross(A,B,2).^2,2)) ./ sqrt(sum(A.^2,2));
end