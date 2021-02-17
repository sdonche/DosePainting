function jsonSARRP(SARRPinput,desired_beams,dose_distri)
%jsonSARRP TODOSummary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_VOI = length(desired_beams);    % Number of VOIs
counter = 1;                        % Counter
fn_VOI = fieldnames(SARRPinput.beam);   % Number of fieldnames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMBER OF BEAMS PER VOI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_beamsVOI = zeros(1,5);
for i = 1:num_VOI
    num_beamsVOI(i) = ...
        length(desired_beams(i).angles_couch0) + ...
        length(desired_beams(i).angles_couchm45) + ...
        length(desired_beams(i).angles_couchm90);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dose = dose_distri ./ num_beamsVOI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERAL (PART 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
json_gen.FileVersion = 5;
json_gen.name = 'Sam';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEAMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(SARRPinput.beam)
    
    angle_tmp = SARRPinput.beam(j).angle;
    couch_tmp = SARRPinput.beam(j).couch_pos;
    
    if couch_tmp == 0
        for K1 = 1:num_VOI
            if any(desired_beams(K1).angles_couch0 == angle_tmp)
                
                    json_gen.studies.beams(counter).adjust_weight = 'Yes';
                    json_gen.studies.beams(counter).collimator = ...
                        strcat('v',num2str(SARRPinput.beam(j).(fn_VOI{K1+3}).('YawDist_X')),'x', ...
                        num2str(SARRPinput.beam(j).(fn_VOI{K1+3}).('YawDist_Y')),' mm');
                    json_gen.studies.beams(counter).collimator_rotation_angle = '0';                % No collimator rotation used
                    json_gen.studies.beams(counter).collimator_x = SARRPinput.beam(j).(fn_VOI{K1+3}).('YawDist_X');
                    json_gen.studies.beams(counter).collimator_y = SARRPinput.beam(j).(fn_VOI{K1+3}).('YawDist_Y');
                    json_gen.studies.beams(counter).couch = couch_tmp;
                    json_gen.studies.beams(counter).gantry = angle_tmp;
                    json_gen.studies.beams(counter).isocenter = strcat('IsoC_',num2str(counter));
                    json_gen.studies.beams(counter).label = strcat('Beam_',num2str(counter));
                    json_gen.studies.beams(counter).ssd = '0';
                    json_gen.studies.beams(counter).time = '0';
                    json_gen.studies.beams(counter).type = 'Beam';
                    json_gen.studies.beams(counter).weight = 100.0;
                    
                    counter = counter + 1;
            end
        end
       
    elseif couch_tmp == -45
        for K2 = 1:num_VOI
            if any(desired_beams(K2).angles_couchm45 == angle_tmp)
                
                    json_gen.studies.beams(counter).adjust_weight = 'Yes';
                    json_gen.studies.beams(counter).collimator = ...
                        strcat('v',num2str(SARRPinput.beam(j).(fn_VOI{K2+3}).('YawDist_X')),'x', ...
                        num2str(SARRPinput.beam(j).(fn_VOI{K2+3}).('YawDist_Y')),' mm');
                    json_gen.studies.beams(counter).collimator_rotation_angle = '0';                % No collimator rotation used
                    json_gen.studies.beams(counter).collimator_x = SARRPinput.beam(j).(fn_VOI{K2+3}).('YawDist_X');
                    json_gen.studies.beams(counter).collimator_y = SARRPinput.beam(j).(fn_VOI{K2+3}).('YawDist_Y');
                    json_gen.studies.beams(counter).couch = couch_tmp;
                    json_gen.studies.beams(counter).gantry = angle_tmp;
                    json_gen.studies.beams(counter).isocenter = strcat('IsoC_',num2str(counter));
                    json_gen.studies.beams(counter).label = strcat('Beam_',num2str(counter));
                    json_gen.studies.beams(counter).ssd = '0';
                    json_gen.studies.beams(counter).time = '0';
                    json_gen.studies.beams(counter).type = 'Beam';
                    json_gen.studies.beams(counter).weight = 100.0;
                    
                    counter = counter + 1;
            end
        end        
        
    elseif couch_tmp == -90 
        for K3 = 1:num_VOI
            if any(desired_beams(K3).angles_couchm90 == angle_tmp)
                
                    json_gen.studies.beams(counter).adjust_weight = 'Yes';
                    json_gen.studies.beams(counter).collimator = ...
                        strcat('v',num2str(SARRPinput.beam(j).(fn_VOI{K3+3}).('YawDist_X')),'x', ...
                        num2str(SARRPinput.beam(j).(fn_VOI{K3+3}).('YawDist_Y')),' mm');
                    json_gen.studies.beams(counter).collimator_rotation_angle = '0';                % No collimator rotation used
                    json_gen.studies.beams(counter).collimator_x = SARRPinput.beam(j).(fn_VOI{K3+3}).('YawDist_X');
                    json_gen.studies.beams(counter).collimator_y = SARRPinput.beam(j).(fn_VOI{K3+3}).('YawDist_Y');
                    json_gen.studies.beams(counter).couch = couch_tmp;
                    json_gen.studies.beams(counter).gantry = angle_tmp;
                    json_gen.studies.beams(counter).isocenter = strcat('IsoC_',num2str(counter));
                    json_gen.studies.beams(counter).label = strcat('Beam_',num2str(counter));
                    json_gen.studies.beams(counter).ssd = '0';
                    json_gen.studies.beams(counter).time = '0';
                    json_gen.studies.beams(counter).type = 'Beam';
                    json_gen.studies.beams(counter).weight = 100.0;
                    
                    counter = counter + 1;
            end
        end              
    else
        
        % TODO ERROR: This couch position is not supported
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
json_gen.studies.contours = [];
json_gen.studies.dose_type = 'Water';
json_gen.studies.heterogene = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Isocenters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter = 1;                % Reset counter
for j = 1:length(SARRPinput.beam)
    
    angle_tmp = SARRPinput.beam(j).angle;
    couch_tmp = SARRPinput.beam(j).couch_pos;
    
    if couch_tmp == 0
        for K1 = 1:num_VOI
            if any(desired_beams(K1).angles_couch0 == angle_tmp)
                
                    json_gen.studies.isocenters(counter).coordinates = SARRPinput.beam(j).(fn_VOI{K1+3}).('YawMid_SARRP');
                    json_gen.studies.isocenters(counter).dose = dose(K1);
                    json_gen.studies.isocenters(counter).name = strcat('IsoC_',num2str(counter));
                    
                    counter = counter + 1;
            end
        end
       
    elseif couch_tmp == -45
        for K2 = 1:num_VOI
            if any(desired_beams(K2).angles_couchm45 == angle_tmp)
                
                    json_gen.studies.isocenters(counter).coordinates = SARRPinput.beam(j).(fn_VOI{K2+3}).('YawMid_SARRP');
                    json_gen.studies.isocenters(counter).dose = dose(K2);
                    json_gen.studies.isocenters(counter).name = strcat('IsoC_',num2str(counter));
                    
                    counter = counter + 1;
            end
        end        
        
    elseif couch_tmp == -90 
        for K3 = 1:num_VOI
            if any(desired_beams(K3).angles_couchm90 == angle_tmp)
                
                    json_gen.studies.isocenters(counter).coordinates = SARRPinput.beam(j).(fn_VOI{K3+3}).('YawMid_SARRP');
                    json_gen.studies.isocenters(counter).dose = dose(K3);
                    json_gen.studies.isocenters(counter).name = strcat('IsoC_',num2str(counter));
                    
                    counter = counter + 1;
            end
        end              
    else
        
        % TODO ERROR: This couch position is not supported
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERAL (PART 3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

json_gen.studies.name = 'Bitcoin_to_the_moon';

%Phantom/false
json_gen.studies.segmentation_thresholds.phantom.false(1).colorName =   'Air';
json_gen.studies.segmentation_thresholds.phantom.false(1).enabled =     true;
json_gen.studies.segmentation_thresholds.phantom.false(1).maximum =     2753;
json_gen.studies.segmentation_thresholds.phantom.false(1).medium =      'Air';
json_gen.studies.segmentation_thresholds.phantom.false(1).minimum =     0;
json_gen.studies.segmentation_thresholds.phantom.false(1).value =       0;

json_gen.studies.segmentation_thresholds.phantom.false(2).colorName =   'Fat';
json_gen.studies.segmentation_thresholds.phantom.false(2).enabled =     true;
json_gen.studies.segmentation_thresholds.phantom.false(2).maximum =     65535;
json_gen.studies.segmentation_thresholds.phantom.false(2).medium =      'Water';
json_gen.studies.segmentation_thresholds.phantom.false(2).minimum =     0;
json_gen.studies.segmentation_thresholds.phantom.false(2).value =       2;

%Phantom/true
json_gen.studies.segmentation_thresholds.phantom.true(1).colorName =    'Air';
json_gen.studies.segmentation_thresholds.phantom.true(1).enabled =      true;
json_gen.studies.segmentation_thresholds.phantom.true(1).maximum =      1298;
json_gen.studies.segmentation_thresholds.phantom.true(1).medium =       'Air';
json_gen.studies.segmentation_thresholds.phantom.true(1).minimum =      0;
json_gen.studies.segmentation_thresholds.phantom.true(1).value =        0;

json_gen.studies.segmentation_thresholds.phantom.true(2).colorName =    'Lung';
json_gen.studies.segmentation_thresholds.phantom.true(2).enabled =      true;
json_gen.studies.segmentation_thresholds.phantom.true(2).maximum =      2753;
json_gen.studies.segmentation_thresholds.phantom.true(2).medium =       'Cork';
json_gen.studies.segmentation_thresholds.phantom.true(2).minimum =      1299;
json_gen.studies.segmentation_thresholds.phantom.true(2).value =        1;

json_gen.studies.segmentation_thresholds.phantom.true(3).colorName =    'Fat';
json_gen.studies.segmentation_thresholds.phantom.true(3).enabled =      true;
json_gen.studies.segmentation_thresholds.phantom.true(3).maximum =      6064;
json_gen.studies.segmentation_thresholds.phantom.true(3).medium =       'Water';
json_gen.studies.segmentation_thresholds.phantom.true(3).minimum =      2754;
json_gen.studies.segmentation_thresholds.phantom.true(3).value =        2;

json_gen.studies.segmentation_thresholds.phantom.true(4).colorName =    'Tissue';
json_gen.studies.segmentation_thresholds.phantom.true(4).enabled =      true;
json_gen.studies.segmentation_thresholds.phantom.true(4).maximum =      9354;
json_gen.studies.segmentation_thresholds.phantom.true(4).medium =       'Graphite';
json_gen.studies.segmentation_thresholds.phantom.true(4).minimum =      6065;
json_gen.studies.segmentation_thresholds.phantom.true(4).value =        3;

json_gen.studies.segmentation_thresholds.phantom.true(5).colorName =    'Bone';
json_gen.studies.segmentation_thresholds.phantom.true(5).enabled =      true;
json_gen.studies.segmentation_thresholds.phantom.true(5).maximum =      65535;
json_gen.studies.segmentation_thresholds.phantom.true(5).medium =       'Aluminum';
json_gen.studies.segmentation_thresholds.phantom.true(5).minimum =      9355;
json_gen.studies.segmentation_thresholds.phantom.true(5).value =        4;

%Tissue/false
json_gen.studies.segmentation_thresholds.tissue.false(1).colorName =    'Air';
json_gen.studies.segmentation_thresholds.tissue.false(1).enabled =      true;
json_gen.studies.segmentation_thresholds.tissue.false(1).maximum =      16100;
json_gen.studies.segmentation_thresholds.tissue.false(1).medium =       'Air';
json_gen.studies.segmentation_thresholds.tissue.false(1).minimum =      0;
json_gen.studies.segmentation_thresholds.tissue.false(1).value =        0;

json_gen.studies.segmentation_thresholds.tissue.false(2).colorName =    'Tissue';
json_gen.studies.segmentation_thresholds.tissue.false(2).enabled =      true;
json_gen.studies.segmentation_thresholds.tissue.false(2).maximum =      65535;
json_gen.studies.segmentation_thresholds.tissue.false(2).medium =       'Tissue';
json_gen.studies.segmentation_thresholds.tissue.false(2).minimum =      16101;
json_gen.studies.segmentation_thresholds.tissue.false(2).value =        3;

%Tissue/true
json_gen.studies.segmentation_thresholds.tissue.true(1).colorName =     'Air';
json_gen.studies.segmentation_thresholds.tissue.true(1).enabled =       true;
json_gen.studies.segmentation_thresholds.tissue.true(1).maximum =       11621;
json_gen.studies.segmentation_thresholds.tissue.true(1).medium =        'Air';
json_gen.studies.segmentation_thresholds.tissue.true(1).minimum =       0;
json_gen.studies.segmentation_thresholds.tissue.true(1).value =         0;

json_gen.studies.segmentation_thresholds.tissue.true(2).colorName =     'Lung';
json_gen.studies.segmentation_thresholds.tissue.true(2).enabled =       true;
json_gen.studies.segmentation_thresholds.tissue.true(2).maximum =       11623;
json_gen.studies.segmentation_thresholds.tissue.true(2).medium =        'Lung';
json_gen.studies.segmentation_thresholds.tissue.true(2).minimum =       11622;
json_gen.studies.segmentation_thresholds.tissue.true(2).value =         1;

json_gen.studies.segmentation_thresholds.tissue.true(3).colorName =     'Fat';
json_gen.studies.segmentation_thresholds.tissue.true(3).enabled =       true;
json_gen.studies.segmentation_thresholds.tissue.true(3).maximum =       11625;
json_gen.studies.segmentation_thresholds.tissue.true(3).medium =        'Fat';
json_gen.studies.segmentation_thresholds.tissue.true(3).minimum =       11624;
json_gen.studies.segmentation_thresholds.tissue.true(3).value =         2;

json_gen.studies.segmentation_thresholds.tissue.true(4).colorName =     'Tissue';
json_gen.studies.segmentation_thresholds.tissue.true(4).enabled =       true;
json_gen.studies.segmentation_thresholds.tissue.true(4).maximum =       12650;
json_gen.studies.segmentation_thresholds.tissue.true(4).medium =        'Fat';
json_gen.studies.segmentation_thresholds.tissue.true(4).minimum =       11626;
json_gen.studies.segmentation_thresholds.tissue.true(4).value =         3;

json_gen.studies.segmentation_thresholds.tissue.true(5).colorName =     'Bone';
json_gen.studies.segmentation_thresholds.tissue.true(5).enabled =       true;
json_gen.studies.segmentation_thresholds.tissue.true(5).maximum =       65535;
json_gen.studies.segmentation_thresholds.tissue.true(5).medium =        'Bone';
json_gen.studies.segmentation_thresholds.tissue.true(5).minimum =       12651;
json_gen.studies.segmentation_thresholds.tissue.true(5).value =         4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segmentation Type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

json_gen.studies.segmentation_type = 'tissue';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADJUSTMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

json_gen.studies.beams = json_gen.studies.beams';
json_gen.studies.isocenters = json_gen.studies.isocenters';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save json file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jsonStr = jsonencode(json_gen);

jsonStr = strrep(jsonStr,'"studies":','"studies":[');   % Insert extra bracket [
jsonStr = strcat(jsonStr(1:end-1),']',jsonStr(end));    % Insert extra bracket ]

fid = fopen('SARRPinput.json', 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);

disp('A json file has been generated with the SARRP input variables in workdirectory.')

% Mogelijke fouten
% - overal spaties tekort


assignin('base','jsontesting',json_gen);

end

