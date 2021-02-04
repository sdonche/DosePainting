function [matrix_normalvect] = generate_normvect_3arcs(angles,couch_pos)
%TODO UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
n = length(angles);
matrix_normalvect = zeros(n,3);
    
for i = 1:n

    angle = angles(i);
    couch = couch_pos(i);
    
    if couch == 0
        
        if angle < 90
            matrix_normalvect(i,2) = tand(angle);
            matrix_normalvect(i,3) = 1;
        elseif angle == 90
            matrix_normalvect(i,2) = 1;
            matrix_normalvect(i,3) = 0;
        elseif angle > 90 && angle <= 120
            matrix_normalvect(i,2) = 1/tand(angle-90);
            matrix_normalvect(i,3) = -1;
        else
            error('Current SARRP setup does not allow angles bigger than 120°.')
        end
        
    elseif couch == -45
        
        if angle < 90
            matrix_normalvect(i,1) = -tand(angle)/tand(couch);
            matrix_normalvect(i,2) = tand(angle);
            matrix_normalvect(i,3) = 1;
        elseif angle == 90
            matrix_normalvect(i,1) = -1/tand(couch);
            matrix_normalvect(i,2) = 1;
            matrix_normalvect(i,3) = 0;
        elseif angle > 90 && angle <= 120
            matrix_normalvect(i,1) = -tand(couch)/tand(angle-90);
            matrix_normalvect(i,2) = -1/tand(angle-90);
            matrix_normalvect(i,3) = -1;
        else
            error('Current SARRP setup does not allow angles bigger than 120°.')
        end
        
    elseif couch == -90
        
        if angle <= 60
            matrix_normalvect(i,1) = tand(angle);
            matrix_normalvect(i,3) = 1;
        else
            error('Current SARRP setup does not allow angles bigger than 60° at couch position -90°.')
        end
        
    else
        error('This function does not allow these couch position(s).')
    end
    
end
