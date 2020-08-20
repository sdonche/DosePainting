%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Q factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% For PTV(69+PET), the resulting dose distribution was compared to the
% intended voxel intensity-based dose pattern; therefore, a quality factor
% (QF) was introduced. In every point p theat is randomly seeded within
% PTV(69+PET) (n points in total), the obtained over-intended dose ratio
% Q(p) = D(p)/(D(I))(p) is calculated. This information can be visualized
% in a Q-volume histogram (QVH), by displaying the partial PTV(69+PET)
% volume for which Q is greater than or equal to each abscis value. Ideally
% such a curve would drop steeply at Q = 1. QF is defined as the mean
% absolute deviation of Q to 1 withing PTV(69+PET):
%
% QF = (1/n)*sum(p)abs(Q(p)-1)