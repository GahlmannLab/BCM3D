function [pG, pS] = CellLocsCaps_CellModeler(radius,lc,sp)
% modified from CellLocsCapsU3
% 7/18/2019 Ji Zhang

%Instead of filling concentric circles in each slice of cell,
%set up fine grid of points, slice through cell at same resolution, and 

%Output
%pG - check what points inside geometry and add to new list
%pS - check what points at the surface of geometry and add to new list


%unit is  nm
r = radius;

[xE,yE,zE]=meshgrid(-5000:sp:5000,-5000:sp:5000,-5000:sp:5000);% this range need chanege according to cells' size
pG = [xE(:),yE(:),zE(:)];
pS = pG;

%End Cap 1 logical
% inSphere1 = sqrt((pG(:,1)+lc/2).^2+pG(:,2).^2+pG(:,3).^2)<=r;
% inSphere1surface = sqrt((pG(:,1)+lc/2).^2+pG(:,2).^2+pG(:,3).^2)>=r-sp;
inSphere1 = (pG(:,1)+lc/2).^2+pG(:,2).^2+pG(:,3).^2 <= r^2;
inSphere1surface = (pG(:,1)+lc/2).^2+pG(:,2).^2+pG(:,3).^2 >= (r-sp)^2;
leftBound = pG(:,1)<=-lc/2;
EC1_inside = inSphere1 & leftBound;
EC1_surface = inSphere1 & leftBound & inSphere1surface;

%End Cap 2 logical
% inSphere2 = sqrt((pG(:,1)-lc/2).^2+pG(:,2).^2+pG(:,3).^2)<=r;
% inSphere2surface = round(sqrt((pG(:,1)-lc/2).^2+pG(:,2).^2+pG(:,3).^2))>=r-sp;
inSphere2 = (pG(:,1)-lc/2).^2+pG(:,2).^2+pG(:,3).^2<=r^2;
inSphere2surface = (pG(:,1)-lc/2).^2+pG(:,2).^2+pG(:,3).^2>=(r-sp)^2;
rightBound = pG(:,1)>=lc/2;
EC2_inside = inSphere2 & rightBound;
EC2_surface = inSphere2 & rightBound & inSphere2surface;

%Cylinder Logical
% inCylinder = sqrt(pG(:,2).^2+pG(:,3).^2)<=r;
% inCylindersurface= sqrt(pG(:,2).^2+pG(:,3).^2)>=r-sp;
inCylinder = pG(:,2).^2+pG(:,3).^2<=r^2;
inCylindersurface= pG(:,2).^2+pG(:,3).^2>=(r-sp)^2;
centerBounds = pG(:,1)>=-lc/2 & pG(:,1)<=lc/2;
C_inside = inCylinder & centerBounds;
C_surface = inCylinder & centerBounds & inCylindersurface;

% output correct points
pG = pG(EC1_inside | EC2_inside | C_inside, 1:3);
pS = pS(EC1_surface | EC2_surface | C_surface, 1:3);


end