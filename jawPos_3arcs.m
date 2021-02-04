function [output] = jawPos_3arcs(angles,couch_pos,centroids,SROW,varargin)
%TODO UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   angles = vector with different angles

%TODO check if number of angles = number of elements in each varargin
%element

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPARATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define variables
num_ang = size(angles,2);                                                   % number of gantry angles TODO niet helemaal correct, eerder aantal beams
num_couch = size(unique(couch_pos),2);                                      % number of couch positions
track_couch = sort(unique(couch_pos),'descend');                            % unique couch positions
num_lvl = size(varargin,2);                                                 % number of dose painting levels
normalvectors = generate_normvect_3arcs(angles,couch_pos);                  % normal vectors
track_voi = [{'VOI50'},{'VOI60'},{'VOI70'},{'VOI80'},{'VOI90'},{'VOI95'}];  % VOIs

% Generate min-max table of projections
% Be careful with adding other terms here, it might ruin the further code
% VOI fieldnames from index 4 atm

for i = 1:num_ang
   
    output.beam(i).angle = angles(i);
    output.beam(i).couch_pos = couch_pos(i);
    output.beam(i).normalvector = normalvectors(i,:);
    
    % VOI50
    output.beam(i).VOI50.Xmin_pix = [];
    output.beam(i).VOI50.Xmax_pix = [];
    output.beam(i).VOI50.Ymin_pix = [];
    output.beam(i).VOI50.Ymax_pix = [];
    output.beam(i).VOI50.Zmin_pix = [];
    output.beam(i).VOI50.Zmax_pix = [];
    output.beam(i).VOI50.BoundingBox_pix.point1 = zeros(1,3);
    output.beam(i).VOI50.BoundingBox_pix.point2 = zeros(1,3);
    output.beam(i).VOI50.BoundingBox_pix.point3 = zeros(1,3);
    output.beam(i).VOI50.BoundingBox_pix.point4 = zeros(1,3);
    output.beam(i).VOI50.BoundingBox_woo.point1 = zeros(1,3);
    output.beam(i).VOI50.BoundingBox_woo.point2 = zeros(1,3);
    output.beam(i).VOI50.BoundingBox_woo.point3 = zeros(1,3);
    output.beam(i).VOI50.BoundingBox_woo.point4 = zeros(1,3);
    output.beam(i).VOI50.BoundingBox_SARRP.point1 = zeros(1,3);
    output.beam(i).VOI50.BoundingBox_SARRP.point2 = zeros(1,3);
    output.beam(i).VOI50.BoundingBox_SARRP.point3 = zeros(1,3);
    output.beam(i).VOI50.BoundingBox_SARRP.point4 = zeros(1,3);
    output.beam(i).VOI50.YawDist_X = [];
    output.beam(i).VOI50.YawDist_Y = [];
    output.beam(i).VOI50.YawMid_woo = zeros(1,3);
    output.beam(i).VOI50.YawMid_SARRP = zeros(1,3);
    output.beam(i).VOI50.Centroid_pix = centroids(1,:);
    output.beam(i).VOI50.Centroid_woo = zeros(1,3);
    output.beam(i).VOI50.Centroid_SARRP = zeros(1,3);
    
    % VOI60
    output.beam(i).VOI60.Xmin_pix = [];
    output.beam(i).VOI60.Xmax_pix = [];
    output.beam(i).VOI60.Ymin_pix = [];
    output.beam(i).VOI60.Ymax_pix = [];
    output.beam(i).VOI60.Zmin_pix = [];
    output.beam(i).VOI60.Zmax_pix = [];
    output.beam(i).VOI60.BoundingBox_pix.point1 = zeros(1,3);
    output.beam(i).VOI60.BoundingBox_pix.point2 = zeros(1,3);
    output.beam(i).VOI60.BoundingBox_pix.point3 = zeros(1,3);
    output.beam(i).VOI60.BoundingBox_pix.point4 = zeros(1,3);
    output.beam(i).VOI60.BoundingBox_woo.point1 = zeros(1,3);
    output.beam(i).VOI60.BoundingBox_woo.point2 = zeros(1,3);
    output.beam(i).VOI60.BoundingBox_woo.point3 = zeros(1,3);
    output.beam(i).VOI60.BoundingBox_woo.point4 = zeros(1,3);
    output.beam(i).VOI60.BoundingBox_SARRP.point1 = zeros(1,3);
    output.beam(i).VOI60.BoundingBox_SARRP.point2 = zeros(1,3);
    output.beam(i).VOI60.BoundingBox_SARRP.point3 = zeros(1,3);
    output.beam(i).VOI60.BoundingBox_SARRP.point4 = zeros(1,3);
    output.beam(i).VOI60.YawDist_X = [];
    output.beam(i).VOI60.YawDist_Y = [];
    output.beam(i).VOI60.YawMid_woo = zeros(1,3);
    output.beam(i).VOI60.YawMid_SARRP = zeros(1,3);
    output.beam(i).VOI60.Centroid_pix = centroids(2,:);
    output.beam(i).VOI60.Centroid_woo = zeros(1,3);
    output.beam(i).VOI60.Centroid_SARRP = zeros(1,3);
    
    % VOI70
    output.beam(i).VOI70.Xmin_pix = [];
    output.beam(i).VOI70.Xmax_pix = [];
    output.beam(i).VOI70.Ymin_pix = [];
    output.beam(i).VOI70.Ymax_pix = [];
    output.beam(i).VOI70.Zmin_pix = [];
    output.beam(i).VOI70.Zmax_pix = [];
    output.beam(i).VOI70.BoundingBox_pix.point1 = zeros(1,3);
    output.beam(i).VOI70.BoundingBox_pix.point2 = zeros(1,3);
    output.beam(i).VOI70.BoundingBox_pix.point3 = zeros(1,3);
    output.beam(i).VOI70.BoundingBox_pix.point4 = zeros(1,3);
    output.beam(i).VOI70.BoundingBox_woo.point1 = zeros(1,3);
    output.beam(i).VOI70.BoundingBox_woo.point2 = zeros(1,3);
    output.beam(i).VOI70.BoundingBox_woo.point3 = zeros(1,3);
    output.beam(i).VOI70.BoundingBox_woo.point4 = zeros(1,3);
    output.beam(i).VOI70.BoundingBox_SARRP.point1 = zeros(1,3);
    output.beam(i).VOI70.BoundingBox_SARRP.point2 = zeros(1,3);
    output.beam(i).VOI70.BoundingBox_SARRP.point3 = zeros(1,3);
    output.beam(i).VOI70.BoundingBox_SARRP.point4 = zeros(1,3);
    output.beam(i).VOI70.YawDist_X = [];
    output.beam(i).VOI70.YawDist_Y = [];
    output.beam(i).VOI70.YawMid_woo = zeros(1,3);
    output.beam(i).VOI70.YawMid_SARRP = zeros(1,3);
    output.beam(i).VOI70.Centroid_pix = centroids(3,:);
    output.beam(i).VOI70.Centroid_woo = zeros(1,3);
    output.beam(i).VOI70.Centroid_SARRP = zeros(1,3);
    
    % VOI80
    output.beam(i).VOI80.Xmin_pix = [];
    output.beam(i).VOI80.Xmax_pix = [];
    output.beam(i).VOI80.Ymin_pix = [];
    output.beam(i).VOI80.Ymax_pix = [];
    output.beam(i).VOI80.Zmin_pix = [];
    output.beam(i).VOI80.Zmax_pix = [];
    output.beam(i).VOI80.BoundingBox_pix.point1 = zeros(1,3);
    output.beam(i).VOI80.BoundingBox_pix.point2 = zeros(1,3);
    output.beam(i).VOI80.BoundingBox_pix.point3 = zeros(1,3);
    output.beam(i).VOI80.BoundingBox_pix.point4 = zeros(1,3);
    output.beam(i).VOI80.BoundingBox_woo.point1 = zeros(1,3);
    output.beam(i).VOI80.BoundingBox_woo.point2 = zeros(1,3);
    output.beam(i).VOI80.BoundingBox_woo.point3 = zeros(1,3);
    output.beam(i).VOI80.BoundingBox_woo.point4 = zeros(1,3);
    output.beam(i).VOI80.BoundingBox_SARRP.point1 = zeros(1,3);
    output.beam(i).VOI80.BoundingBox_SARRP.point2 = zeros(1,3);
    output.beam(i).VOI80.BoundingBox_SARRP.point3 = zeros(1,3);
    output.beam(i).VOI80.BoundingBox_SARRP.point4 = zeros(1,3);
    output.beam(i).VOI80.YawDist_X = [];
    output.beam(i).VOI80.YawDist_Y = [];
    output.beam(i).VOI80.YawMid_woo = zeros(1,3);
    output.beam(i).VOI80.YawMid_SARRP = zeros(1,3);
    output.beam(i).VOI80.Centroid_pix = centroids(4,:);
    output.beam(i).VOI80.Centroid_woo = zeros(1,3);
    output.beam(i).VOI80.Centroid_SARRP = zeros(1,3);
    
    % VOI90
    output.beam(i).VOI90.Xmin_pix = [];
    output.beam(i).VOI90.Xmax_pix = [];
    output.beam(i).VOI90.Ymin_pix = [];
    output.beam(i).VOI90.Ymax_pix = [];
    output.beam(i).VOI90.Zmin_pix = [];
    output.beam(i).VOI90.Zmax_pix = [];
    output.beam(i).VOI90.BoundingBox_pix.point1 = zeros(1,3);
    output.beam(i).VOI90.BoundingBox_pix.point2 = zeros(1,3);
    output.beam(i).VOI90.BoundingBox_pix.point3 = zeros(1,3);
    output.beam(i).VOI90.BoundingBox_pix.point4 = zeros(1,3);
    output.beam(i).VOI90.BoundingBox_woo.point1 = zeros(1,3);
    output.beam(i).VOI90.BoundingBox_woo.point2 = zeros(1,3);
    output.beam(i).VOI90.BoundingBox_woo.point3 = zeros(1,3);
    output.beam(i).VOI90.BoundingBox_woo.point4 = zeros(1,3);
    output.beam(i).VOI90.BoundingBox_SARRP.point1 = zeros(1,3);
    output.beam(i).VOI90.BoundingBox_SARRP.point2 = zeros(1,3);
    output.beam(i).VOI90.BoundingBox_SARRP.point3 = zeros(1,3);
    output.beam(i).VOI90.BoundingBox_SARRP.point4 = zeros(1,3);
    output.beam(i).VOI90.YawDist_X = [];
    output.beam(i).VOI90.YawDist_Y = [];
    output.beam(i).VOI90.YawMid_woo = zeros(1,3);
    output.beam(i).VOI90.YawMid_SARRP = zeros(1,3);
    output.beam(i).VOI90.Centroid_pix = centroids(5,:);
    output.beam(i).VOI90.Centroid_woo = zeros(1,3);
    output.beam(i).VOI90.Centroid_SARRP = zeros(1,3);
    
    % VOI95
    output.beam(i).VOI95.Xmin_pix = [];
    output.beam(i).VOI95.Xmax_pix = [];
    output.beam(i).VOI95.Ymin_pix = [];
    output.beam(i).VOI95.Ymax_pix = [];
    output.beam(i).VOI95.Zmin_pix = [];
    output.beam(i).VOI95.Zmax_pix = [];
    output.beam(i).VOI95.BoundingBox_pix.point1 = zeros(1,3);
    output.beam(i).VOI95.BoundingBox_pix.point2 = zeros(1,3);
    output.beam(i).VOI95.BoundingBox_pix.point3 = zeros(1,3);
    output.beam(i).VOI95.BoundingBox_pix.point4 = zeros(1,3);
    output.beam(i).VOI95.BoundingBox_woo.point1 = zeros(1,3);
    output.beam(i).VOI95.BoundingBox_woo.point2 = zeros(1,3);
    output.beam(i).VOI95.BoundingBox_woo.point3 = zeros(1,3);
    output.beam(i).VOI95.BoundingBox_woo.point4 = zeros(1,3); 
    output.beam(i).VOI95.BoundingBox_SARRP.point1 = zeros(1,3);
    output.beam(i).VOI95.BoundingBox_SARRP.point2 = zeros(1,3);
    output.beam(i).VOI95.BoundingBox_SARRP.point3 = zeros(1,3);
    output.beam(i).VOI95.BoundingBox_SARRP.point4 = zeros(1,3);
    output.beam(i).VOI95.YawDist_X = [];
    output.beam(i).VOI95.YawDist_Y = [];
    output.beam(i).VOI95.YawMid_woo = zeros(1,3);
    output.beam(i).VOI95.YawMid_SARRP = zeros(1,3);
    output.beam(i).VOI95.Centroid_pix = centroids(6,:);
    output.beam(i).VOI95.Centroid_woo = zeros(1,3);
    output.beam(i).VOI95.Centroid_SARRP = zeros(1,3);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXTRACT MIN/MAX VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TODO make it check whether or not it's really VOI50, not just k=1, to
%avoid mistakes in future

for j = 1:num_ang
    for k = 1:num_lvl
   
    min_tmp = min(varargin{k}{j});
    max_tmp = max(varargin{k}{j});
    
        if k == 1       % VOI50 TODO make this a check for 'VOI50' instead of equal to 1
            output.beam(j).VOI50.Xmax_pix = max_tmp(1,1);
            output.beam(j).VOI50.Xmin_pix = min_tmp(1,1);
            output.beam(j).VOI50.Ymax_pix = max_tmp(1,2);
            output.beam(j).VOI50.Ymin_pix = min_tmp(1,2);
            output.beam(j).VOI50.Zmax_pix = max_tmp(1,3);
            output.beam(j).VOI50.Zmin_pix = min_tmp(1,3);
    
        elseif k == 2   % VOI60
            output.beam(j).VOI60.Xmax_pix = max_tmp(1,1);
            output.beam(j).VOI60.Xmin_pix = min_tmp(1,1);
            output.beam(j).VOI60.Ymax_pix = max_tmp(1,2);
            output.beam(j).VOI60.Ymin_pix = min_tmp(1,2);
            output.beam(j).VOI60.Zmax_pix = max_tmp(1,3);
            output.beam(j).VOI60.Zmin_pix = min_tmp(1,3);
            
        elseif k == 3   % VOI70
            output.beam(j).VOI70.Xmax_pix = max_tmp(1,1);
            output.beam(j).VOI70.Xmin_pix = min_tmp(1,1);
            output.beam(j).VOI70.Ymax_pix = max_tmp(1,2);
            output.beam(j).VOI70.Ymin_pix = min_tmp(1,2);
            output.beam(j).VOI70.Zmax_pix = max_tmp(1,3);
            output.beam(j).VOI70.Zmin_pix = min_tmp(1,3);
            
        elseif k == 4   % VOI80
            output.beam(j).VOI80.Xmax_pix = max_tmp(1,1);
            output.beam(j).VOI80.Xmin_pix = min_tmp(1,1);
            output.beam(j).VOI80.Ymax_pix = max_tmp(1,2);
            output.beam(j).VOI80.Ymin_pix = min_tmp(1,2);
            output.beam(j).VOI80.Zmax_pix = max_tmp(1,3);
            output.beam(j).VOI80.Zmin_pix = min_tmp(1,3);
            
        elseif k == 5   % VOI90
            output.beam(j).VOI90.Xmax_pix = max_tmp(1,1);
            output.beam(j).VOI90.Xmin_pix = min_tmp(1,1);
            output.beam(j).VOI90.Ymax_pix = max_tmp(1,2);
            output.beam(j).VOI90.Ymin_pix = min_tmp(1,2);
            output.beam(j).VOI90.Zmax_pix = max_tmp(1,3);
            output.beam(j).VOI90.Zmin_pix = min_tmp(1,3);
            
        elseif k == 6   % VOI95
            output.beam(j).VOI95.Xmax_pix = max_tmp(1,1);
            output.beam(j).VOI95.Xmin_pix = min_tmp(1,1);
            output.beam(j).VOI95.Ymax_pix = max_tmp(1,2);
            output.beam(j).VOI95.Ymin_pix = min_tmp(1,2);
            output.beam(j).VOI95.Zmax_pix = max_tmp(1,3);
            output.beam(j).VOI95.Zmin_pix = min_tmp(1,3);
            
        else
            error('Function error. Investigate!')
        
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BOUNDING BOX (PIXEL COORDINATES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn_VOI = fieldnames(output.beam(1));

for l = 1:num_ang
    for m = 4:length(fn_VOI)            %VOI fieldnames start at 4, careful with adding extra fields to beams!!
    angle_tmp = output.beam(l).angle;
    couch_tmp = output.beam(l).couch_pos;
    %VOI_tmp = output.beam(l).fn_VOI(m);
        
        if angle_tmp == 0                       % Gantry angle 0°
            if couch_tmp == 0                   % Couch angle 0° (default position)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                % Point 1
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point1') = [...
                    output.beam(l).(fn_VOI{m}).('Xmin_pix'),...   % X-value
                    output.beam(l).(fn_VOI{m}).('Ymin_pix'),...   % Y-value
                    output.beam(l).(fn_VOI{m}).('Zmin_pix')];     % Z-value; Zmin = Zmax for gantry angle 0°

                % Point 2
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point2') = [...
                    output.beam(l).(fn_VOI{m}).('Xmin_pix'),...   % X-value
                    output.beam(l).(fn_VOI{m}).('Ymax_pix'),...   % Y-value
                    output.beam(l).(fn_VOI{m}).('Zmin_pix')];     % Z-value; Zmin = Zmax for gantry angle 0°
            
                % Point 3
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point3') = [...
                    output.beam(l).(fn_VOI{m}).('Xmax_pix'),...   % X-value
                    output.beam(l).(fn_VOI{m}).('Ymin_pix'),...   % Y-value
                    output.beam(l).(fn_VOI{m}).('Zmin_pix')];     % Z-value; Zmin = Zmax for gantry angle 0°
            
                % Point 4
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point4') = [...
                    output.beam(l).(fn_VOI{m}).('Xmax_pix'),...   % X-value
                    output.beam(l).(fn_VOI{m}).('Ymax_pix'),...   % Y-value
                    output.beam(l).(fn_VOI{m}).('Zmin_pix')];     % Z-value; Zmin = Zmax for gantry angle 0°
                
            elseif couch_tmp == -90             % Couch angle -90°
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Point 1
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point1') = [...
                    output.beam(l).(fn_VOI{m}).('Xmin_pix'),...   % X-value
                    output.beam(l).(fn_VOI{m}).('Ymin_pix'),...   % Y-value
                    output.beam(l).(fn_VOI{m}).('Zmin_pix')];     % Z-value; Zmin = Zmax for gantry angle 0°

                % Point 2
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point3') = [...
                    output.beam(l).(fn_VOI{m}).('Xmin_pix'),...   % X-value
                    output.beam(l).(fn_VOI{m}).('Ymax_pix'),...   % Y-value
                    output.beam(l).(fn_VOI{m}).('Zmin_pix')];     % Z-value; Zmin = Zmax for gantry angle 0°
            
                % Point 3
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point2') = [...
                    output.beam(l).(fn_VOI{m}).('Xmax_pix'),...   % X-value
                    output.beam(l).(fn_VOI{m}).('Ymin_pix'),...   % Y-value
                    output.beam(l).(fn_VOI{m}).('Zmin_pix')];     % Z-value; Zmin = Zmax for gantry angle 0°
            
                % Point 4
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point4') = [...
                    output.beam(l).(fn_VOI{m}).('Xmax_pix'),...   % X-value
                    output.beam(l).(fn_VOI{m}).('Ymax_pix'),...   % Y-value
                    output.beam(l).(fn_VOI{m}).('Zmin_pix')];     % Z-value; Zmin = Zmax for gantry angle 0°                
                
            else                                % Couch angle < 0° & > -90° (e.g. -45°)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %TODO these need a collimator angle!!!! 
                % E.g. couch = -45, then a collimator angle of 45 is needed
                % distances are the same as couch position 0°
                
                % Point 1
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point1') = [...
                    output.beam(l).(fn_VOI{m}).('Xmin_pix'),...   % X-value
                    output.beam(l).(fn_VOI{m}).('Ymin_pix'),...   % Y-value
                    output.beam(l).(fn_VOI{m}).('Zmin_pix')];     % Z-value; Zmin = Zmax for gantry angle 0°

                % Point 2
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point2') = [...
                    output.beam(l).(fn_VOI{m}).('Xmin_pix'),...   % X-value
                    output.beam(l).(fn_VOI{m}).('Ymax_pix'),...   % Y-value
                    output.beam(l).(fn_VOI{m}).('Zmin_pix')];     % Z-value; Zmin = Zmax for gantry angle 0°
            
                % Point 3
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point3') = [...
                    output.beam(l).(fn_VOI{m}).('Xmax_pix'),...   % X-value
                    output.beam(l).(fn_VOI{m}).('Ymin_pix'),...   % Y-value
                    output.beam(l).(fn_VOI{m}).('Zmin_pix')];     % Z-value; Zmin = Zmax for gantry angle 0°
            
                % Point 4
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point4') = [...
                    output.beam(l).(fn_VOI{m}).('Xmax_pix'),...   % X-value
                    output.beam(l).(fn_VOI{m}).('Ymax_pix'),...   % Y-value
                    output.beam(l).(fn_VOI{m}).('Zmin_pix')];     % Z-value; Zmin = Zmax for gantry angle 0°                
                
            end
        elseif angle_tmp == 90                  % Gantry angle 90°            
            error('Implement exception 2 (for gantry angles at 90°)!!!')
            %TODO
            
            %     elseif angles(k) == 90  % Exception 2
            %         for l = 1:num_lvl
            % 
            %             % Calculate corner points for exception 1
            %             CornerPoints = [min_max_values.(fn{k}){l,1} min_max_values.(fn{k}){l,3} min_max_values.(fn{k}){l,5};...
            %                 min_max_values.(fn{k}){l,1} min_max_values.(fn{k}){l,3} min_max_values.(fn{k}){l,6};...
            %                 min_max_values.(fn{k}){l,2} min_max_values.(fn{k}){l,3} min_max_values.(fn{k}){l,5};...
            %                 min_max_values.(fn{k}){l,2} min_max_values.(fn{k}){l,3} min_max_values.(fn{k}){l,6}];
            %             
            %            % Ymin = Ymax; min_max_values.(fn{k}){l,3} = min_max_values.(fn{k}){l,4}
            %             
            %            Bound_Pix{:,l+1}(tracker(k):tracker(k)+3,:) = CornerPoints;
            %         end   
            if couch_tmp == 0                   % Couch angle 0° (default position) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            elseif couch_tmp == -90             % Couch angle -90°
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            else                                % Couch angle < 0° & > -90° (e.g. -45°)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end
        else                                    % Gantry angle > 0° & < 90° or > 90° & =< 120°
            % INFO: use X-& Y-values or Y-& Z-values from min/max calculation AND
            % f-values to determine Z- or X-values from the bounding box
            if couch_tmp == 0                   % Couch angle 0° (default position)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Calculate f-values
                f1 = output.beam(l).normalvector(1);
                f2 = output.beam(l).normalvector(2);
                f3 = output.beam(l).normalvector(3);
                f4 = -output.beam(l).normalvector(1)*output.beam(l).(fn_VOI{m}).('Centroid_pix')(1)...
                    -output.beam(l).normalvector(2)*output.beam(l).(fn_VOI{m}).('Centroid_pix')(2)...
                    -output.beam(l).normalvector(3)*output.beam(l).(fn_VOI{m}).('Centroid_pix')(3);
            
                % Extract X and Y values
                BoundingPts_pix = [output.beam(l).(fn_VOI{m}).('Xmin_pix') output.beam(l).(fn_VOI{m}).('Ymin_pix') 0;...
                    output.beam(l).(fn_VOI{m}).('Xmin_pix') output.beam(l).(fn_VOI{m}).('Ymax_pix') 0;...
                    output.beam(l).(fn_VOI{m}).('Xmax_pix') output.beam(l).(fn_VOI{m}).('Ymin_pix') 0;...
                    output.beam(l).(fn_VOI{m}).('Xmax_pix') output.beam(l).(fn_VOI{m}).('Ymax_pix') 0];
                
                % Calculate Z values
                Z = -(f1.*BoundingPts_pix(:,1) + f2.*BoundingPts_pix(:,2) + f4)/f3;
                BoundingPts_pix(:,3) = Z;
            
                % Store in output
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point1') = BoundingPts_pix(1,:);   % Point 1
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point2') = BoundingPts_pix(2,:);   % Point 2
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point3') = BoundingPts_pix(3,:);   % Point 3
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point4') = BoundingPts_pix(4,:);   % Point 4 
                
            elseif couch_tmp == -90             % Couch angle -90°
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                % veranderd: bij store in output
                % punt 2 en 3 omgewisseld
                % Calculate f-values
                f1 = output.beam(l).normalvector(1);
                f2 = output.beam(l).normalvector(2);
                f3 = output.beam(l).normalvector(3);
                f4 = -output.beam(l).normalvector(1)*output.beam(l).(fn_VOI{m}).('Centroid_pix')(1)...
                    -output.beam(l).normalvector(2)*output.beam(l).(fn_VOI{m}).('Centroid_pix')(2)...
                    -output.beam(l).normalvector(3)*output.beam(l).(fn_VOI{m}).('Centroid_pix')(3);
            
                % Extract X and Y values
                BoundingPts_pix = [output.beam(l).(fn_VOI{m}).('Xmin_pix') output.beam(l).(fn_VOI{m}).('Ymin_pix') 0;...
                    output.beam(l).(fn_VOI{m}).('Xmin_pix') output.beam(l).(fn_VOI{m}).('Ymax_pix') 0;...
                    output.beam(l).(fn_VOI{m}).('Xmax_pix') output.beam(l).(fn_VOI{m}).('Ymin_pix') 0;...
                    output.beam(l).(fn_VOI{m}).('Xmax_pix') output.beam(l).(fn_VOI{m}).('Ymax_pix') 0];
                
                % Calculate Z values
                Z = -(f1*BoundingPts_pix(:,1) + f2*BoundingPts_pix(:,2) + f4)/f3;
                BoundingPts_pix(:,3) = Z;
            
                % Store in output
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point1') = BoundingPts_pix(1,:);   % Point 1
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point3') = BoundingPts_pix(2,:);   % Point 2
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point2') = BoundingPts_pix(3,:);   % Point 3
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point4') = BoundingPts_pix(4,:);   % Point 4 
                
            else                                % Couch angle < 0° & > -90° (e.g. -45°)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                % Veranderd: BoundingPts_pix helemaal veranderd: juiste
                % volgorde?
                % Ook calculate X values ipv Z values
                % Calculate f-values
                f1 = output.beam(l).normalvector(1);
                f2 = output.beam(l).normalvector(2);
                f3 = output.beam(l).normalvector(3);
                f4 = -output.beam(l).normalvector(1)*output.beam(l).(fn_VOI{m}).('Centroid_pix')(1)...
                    -output.beam(l).normalvector(2)*output.beam(l).(fn_VOI{m}).('Centroid_pix')(2)...
                    -output.beam(l).normalvector(3)*output.beam(l).(fn_VOI{m}).('Centroid_pix')(3);
            
                % Extract Y-and Z-values
                BoundingPts_pix = [0 output.beam(l).(fn_VOI{m}).('Ymin_pix') output.beam(l).(fn_VOI{m}).('Zmin_pix');...
                    0 output.beam(l).(fn_VOI{m}).('Ymin_pix') output.beam(l).(fn_VOI{m}).('Zmax_pix');...
                    0 output.beam(l).(fn_VOI{m}).('Ymax_pix') output.beam(l).(fn_VOI{m}).('Zmin_pix');...
                    0 output.beam(l).(fn_VOI{m}).('Ymax_pix') output.beam(l).(fn_VOI{m}).('Zmax_pix')];
                
                % Calculate X values
                X = -(f2*BoundingPts_pix(:,2) + f3*BoundingPts_pix(:,3) + f4)/f1;
                BoundingPts_pix(:,1) = X;
            
                % Store in output
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point1') = BoundingPts_pix(1,:);   % Point 1
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point3') = BoundingPts_pix(2,:);   % Point 2
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point2') = BoundingPts_pix(3,:);   % Point 3
                output.beam(l).(fn_VOI{m}).('BoundingBox_pix').('point4') = BoundingPts_pix(4,:);   % Point 4                 
            end
    end
end


%TODO are all corner points correct? CHECK -45° couch position!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WORLD COORDINATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% World coordinates are NOT SARRP coordinates!!!
% Perform the transformation below:
% X_woo = Z_SARRP
% Y_woo = - Y_SARRP
% Z_woo = X_SARRP

coor_bucket = zeros(4,4);
coor_bucket(:,4) = 1;
cent_bucket = zeros(1,4);
cent_bucket(1,4) = 1;

% Calculate bounding box points + centroids
for n = 1:num_ang
    for p = 4:length(fn_VOI)            %VOI fieldnames start at 4, careful with adding extra fields to beams!!
        
        % Bounding Box Points
        coor_bucket(1,1:3) = output.beam(n).(fn_VOI{p}).('BoundingBox_pix').('point1');         % Point 1
        coor_bucket(2,1:3) = output.beam(n).(fn_VOI{p}).('BoundingBox_pix').('point2');         % Point 2
        coor_bucket(3,1:3) = output.beam(n).(fn_VOI{p}).('BoundingBox_pix').('point3');         % Point 3
        coor_bucket(4,1:3) = output.beam(n).(fn_VOI{p}).('BoundingBox_pix').('point4');         % Point 4
        
        % Centroid point
        cent_bucket(1,1:3) = output.beam(n).(fn_VOI{p}).('Centroid_pix');
        
        % Pixel to World coordinates
        coor_bucket_tmp = SROW * coor_bucket';
        cent_bucket_tmp = SROW * cent_bucket';
        
        % Store into output variable
        output.beam(n).(fn_VOI{p}).('BoundingBox_woo').('point1') = coor_bucket_tmp(1:3,1)';    % Point 1 
        output.beam(n).(fn_VOI{p}).('BoundingBox_woo').('point2') = coor_bucket_tmp(1:3,2)';    % Point 2
        output.beam(n).(fn_VOI{p}).('BoundingBox_woo').('point3') = coor_bucket_tmp(1:3,3)';    % Point 3
        output.beam(n).(fn_VOI{p}).('BoundingBox_woo').('point4') = coor_bucket_tmp(1:3,4)';    % Point 4
        output.beam(n).(fn_VOI{p}).('Centroid_woo') = cent_bucket_tmp(1:3,1)';                  % Centroid
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE SARRP SETTINGS (WORLD COORDINATES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO moeten de jaws worden gekanteld bij de verschillende couch
% posities??????

% Yaw distances + Yaw mid
for q = 1:num_ang
    for r = 4:length(fn_VOI)            %VOI fieldnames start at 4, careful with adding extra fields to beams!!
        if output.beam(q).couch_pos == -45
            % Distance calculation for Yaw_dist_1, Yaw_dist_2 and Yaw_mid
            % have to be adjusted for couch_pos -45
            Yaw_dist_1 = point2line(output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point3'),...
                output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point1'),...
                output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point2'));
            Yaw_dist_2 = sqrt(sum((output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point1') - ...
                output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point2')).^2 ));
            % Yaw mid
            Yaw_mid = (output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point1') + ...
                output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point4'))/2;
            
            % Store into output
            output.beam(q).(fn_VOI{r}).('YawDist_X') = round(Yaw_dist_2,2);     % Yaw distance X
            output.beam(q).(fn_VOI{r}).('YawDist_Y') = round(Yaw_dist_1,2);     % Yaw distance Y
            output.beam(q).(fn_VOI{r}).('YawMid_woo') = round(Yaw_mid,2)';           % Yaw mid
            
        else
            % Distance point 1 - point 2; sqrt(sum((point1 - point2) .^ 2))
            Yaw_dist_1 = sqrt(sum((output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point1') - ...
                output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point2')).^2 ));

            % Distance point 1 - point 3, distance between point 2 - is diagonal
            Yaw_dist_2 = sqrt(sum((output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point1') - ...
                output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point3')).^2 ));

            % Yaw mid
            % The centroid of any n points is simply the average of the vectors defining them.
            Yaw_mid = (output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point1')' + ...
                output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point2')' + ...
                output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point3')' + ...
                output.beam(q).(fn_VOI{r}).('BoundingBox_woo').('point4')')/4;

            % Store into output
            output.beam(q).(fn_VOI{r}).('YawDist_X') = round(Yaw_dist_2,2);     % Yaw distance X
            output.beam(q).(fn_VOI{r}).('YawDist_Y') = round(Yaw_dist_1,2);     % Yaw distance Y
            output.beam(q).(fn_VOI{r}).('YawMid_woo') = round(Yaw_mid,2)';           % Yaw mid
            
        end  
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SARRP COORDINATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Translate world coordinates to SARRP coordinates
% X_woo = Z_SARRP
% Y_woo = - Y_SARRP
% Z_woo = X_SARRP
for s = 1:num_ang
    for t = 4:length(fn_VOI)            %VOI fieldnames start at 4, careful with adding extra fields to beams!!
        % Bounding Box Points
        A1 = output.beam(s).(fn_VOI{t}).('BoundingBox_woo').('point1');
        A2 = output.beam(s).(fn_VOI{t}).('BoundingBox_woo').('point2');
        A3 = output.beam(s).(fn_VOI{t}).('BoundingBox_woo').('point3');
        A4 = output.beam(s).(fn_VOI{t}).('BoundingBox_woo').('point4');
        
        % Centroid
        B = output.beam(s).(fn_VOI{t}).('Centroid_woo');

        % Yaw mid
        C = output.beam(s).(fn_VOI{t}).('YawMid_woo');
        
        % Store back into output
        output.beam(s).(fn_VOI{t}).('BoundingBox_SARRP').('point1') = [A1(3) -A1(2) A1(1)];
        output.beam(s).(fn_VOI{t}).('BoundingBox_SARRP').('point2') = [A2(3) -A2(2) A2(1)];
        output.beam(s).(fn_VOI{t}).('BoundingBox_SARRP').('point3') = [A3(3) -A3(2) A3(1)];
        output.beam(s).(fn_VOI{t}).('BoundingBox_SARRP').('point4') = [A4(3) -A4(2) A4(1)];
        output.beam(s).(fn_VOI{t}).('Centroid_SARRP') = [B(3) -B(2) B(1)];
        output.beam(s).(fn_VOI{t}).('YawMid_SARRP') = [C(3) -C(2) C(1)];
        
    end
end
end