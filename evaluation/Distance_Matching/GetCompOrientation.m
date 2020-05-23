%% tools
function direction = GetCompOrientation(currentSeg)
% This function is to find the direction of current component
Sub = ones(size(currentSeg));
Sub = mean(currentSeg,1).*Sub;
ShiftSeg = currentSeg - Sub;

CovM = cov(ShiftSeg);
[V,D] = eig(CovM);
x_v1= V(1,3);
y_v1= V(2,3);
z_v1= V(3,3);
direction = [x_v1,y_v1,z_v1];
end