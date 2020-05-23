function [pG, pS] = CellLocsCaps_CellModeler(radius,lc,sp)
%Set up fine grid of points to represent a cell
%Input
%radius - radius of a cell
%length - length of a cell
%sp - grid unit size

%Output
%pG - check what points inside geometry and add to new list
%pS - check what points at the surface of geometry and add to new list


r = radius;

[xE,yE,zE]=meshgrid(-4000:sp:4000,-4000:sp:4000,-4000:sp:4000);% this range need chanege according to cells' size
pG = [xE(:),yE(:),zE(:)];
pS = pG;

%End Cap 1 logical
inSphere1 = sqrt((pG(:,1)+lc/2).^2+pG(:,2).^2+pG(:,3).^2)<=r;
inSphere1surface = sqrt((pG(:,1)+lc/2).^2+pG(:,2).^2+pG(:,3).^2)>=r-sp;
leftBound = pG(:,1)<=-lc/2;
EC1_inside = inSphere1 & leftBound;
EC1_surface = inSphere1 & leftBound & inSphere1surface;

%End Cap 2 logical
inSphere2 = sqrt((pG(:,1)-lc/2).^2+pG(:,2).^2+pG(:,3).^2)<=r;
inSphere2surface = round(sqrt((pG(:,1)-lc/2).^2+pG(:,2).^2+pG(:,3).^2))>=r-sp;
rightBound = pG(:,1)>=lc/2;
EC2_inside = inSphere2 & rightBound;
EC2_surface = inSphere2 & rightBound & inSphere2surface;

%Cylinder Logical
inCylinder = sqrt(pG(:,2).^2+pG(:,3).^2)<=r;
inCylindersurface= sqrt(pG(:,2).^2+pG(:,3).^2)>=r-sp;
centerBounds = pG(:,1)>=-lc/2 & pG(:,1)<=lc/2;
C_inside = inCylinder & centerBounds;
C_surface = inCylinder & centerBounds & inCylindersurface;

% output correct points
pG = pG(EC1_inside | EC2_inside | C_inside, 1:3);
pS = pS(EC1_surface | EC2_surface | C_surface, 1:3);


end