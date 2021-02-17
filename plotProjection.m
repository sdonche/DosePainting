function HAxes = plotProjection(normvector,centroid,original,projection)
%TODO UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%TODO also write equation plane here etc

% Define parameters
f1 = normvector(1);
f2 = normvector(2);
f3 = normvector(3);
f4 = -normvector(1)*centroid(1)-normvector(2)*centroid(2)-normvector(3)*centroid(3);

% Initiate figure
figure('Name','','NumberTitle','off')
HAxes = axes('NextPlot', 'add');
xlabel('HEAD <-- X-as --> TAIL')
ylabel('RIGHT <-- Y-as --> LEFT')
zlabel('DORSAL <-- Z-as --> VENTRAL')
xlim([0 500]) 
ylim([0 500])
zlim([0 500])

% Plot plane
if isempty(normvector) == 0
    [x, y] = meshgrid(1:5:500);                                     % Generate x and y data
    z = -1./f3*(f1.*x + f2.*y + f4);                                % Solve for z data
    surf(x,y,z,'EdgeColor','k','FaceColor','none','Parent',HAxes)   %Plot the surface

    % Incident beam
    t = -1000:1:1000;                                               % Generate "points"
    Xl = centroid(1) + t * normvector(1);
    Yl = centroid(2) + t * normvector(2);
    Zl = centroid(3) + t * normvector(3);
    plot3(Xl,Yl,Zl,'.m','Parent',HAxes)                             % Line through isocenter, orthogonal to projection plane
end

% Plot tumour
if isempty(original) == 0
    scatter3(original(:,1),original(:,2),original(:,3),'ob','Parent',HAxes)
end

%Plot tumour projection
if isempty(projection) == 0
    scatter3(projection(:,1),projection(:,2),projection(:,3),'or','Parent',HAxes)
end
end

