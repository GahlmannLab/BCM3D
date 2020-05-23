function [numCells,rawLocs,ground,surface_select_points,whole_select_points,...
    ground_label] = cellVolume_cellmodeller(CellModeller_data,voxelSize,fov, Distance_ratio, n)

%Pos/orient cells randomly in voxel grid
%Uses larger voxel subgrid to setup regions 
%typical sp = 30 nm ??
%
%colTol sets amount of overlap between cells in #voxels set to 0 don't
%allow overlay
% randomly choose 500-1000 points from a cell 

%g3DLocs-data for each cell: center point, theta, phi, r, lc(temporarily
%not available)

%numCells-number of cells
%MeanPoints - mean number of points selected from each cell
%rawLocs-coordinate for every point
%ground-the binary image for cells
%select_points-points select from each cells' surface used for convolution
%ground_label: cell image with annotations

%g3DLocs=cell(1,5);
xLocs=[];
yLocs=[];
zLocs=[];
numCells = size(CellModeller_data,1);

% translation matrix
COM = Distance_ratio.*CellModeller_data(:, 4:6) + 0.5 * [fov(1),fov(2),1000];%0.5 * [fov(1),fov(2),200];
% ehance by distance _ratio

% cell orientation
cell_direction_vectors = CellModeller_data(:, 7:9);
% cell radius and length
r = CellModeller_data(:, 2);
lc = CellModeller_data(:, 3);

sp=30;%??

%Creating voxel grid w/ resolution res
%res=100;
xRange = 0:voxelSize:fov(1);
yRange = 0:voxelSize:fov(2);
zRange = 0:voxelSize:fov(3);

%preallocate
ground = zeros(length(xRange),length(yRange),length(zRange));
ground_label = ground;
surface_select_points = zeros(length(xRange),length(yRange),length(zRange));
whole_select_points = surface_select_points;

%number of points chosen from each cell for convolution
num_points=zeros(numCells, 1);
for i=1:numCells
    num_points(i)=round(n + rand(1)*500);% change to n - n+500;% old number of points for each cell is from 500-1000
end

%% loop through the CellModeller result to generate all cells
for i = 1:numCells
%     alpha(1) = rand(1)*pi;
%     theta(1) = rand(1)*pi;
%     psi(1) = rand(1)*pi;
    dx = COM(i,1);
    dy = COM(i,2);
    dz = COM(i,3);
    
    %% generate points for the cell
    [All_points,surface_points] = CellLocsCaps_CellModeler(r(i),lc(i),sp);
    
    original_unit_vector = [1 0 0];
    final_unit_vector = cell_direction_vectors(i, :);
    rotated_cell_points = rotCellWithVector(All_points, original_unit_vector, final_unit_vector);
    
    cellXs = rotated_cell_points(:,1);
    cellYs = rotated_cell_points(:,2);
    cellZs = rotated_cell_points(:,3);
    
    % move the cell
    cellXs = cellXs + dx;
    cellYs = cellYs + dy;
    cellZs = cellZs + dz;
    
    % store cell location
    xLocs = [xLocs; cellXs];
    yLocs = [yLocs; cellYs];
    zLocs = [zLocs; cellZs];
    
    %% generate surface points for the cell
    rotated_cell_surface_points = rotCellWithVector(surface_points, original_unit_vector, final_unit_vector);
    
    cellXs_surface = rotated_cell_surface_points(:,1);
    cellYs_surface = rotated_cell_surface_points(:,2);
    cellZs_surface = rotated_cell_surface_points(:,3);
    
    % move the cell surface
    cellXs_surface = cellXs_surface + dx;
    cellYs_surface = cellYs_surface + dy;
    cellZs_surface = cellZs_surface + dz;
    
    %% Converting continuous cell locs to voxel coordinates and filling ground
    %Preallocate
    xC = zeros(length(cellXs),1);
    yC = xC;
    zC = xC;
    for n=1:length(cellXs)
        xC(n) = ceil(cellXs(n)/voxelSize);
        yC(n) = length(yRange)-ceil(cellYs(n)/voxelSize); % 
        zC(n) = ceil(cellZs(n)/voxelSize);
        try
            ground(xC(n),yC(n),zC(n))= 1;
            ground_label(xC(n),yC(n),zC(n))= i;
        catch
            warning('Problem of indexing 1.');
            disp(zC(n));
            ground(xC(n),yC(n),1)= 1;
            ground_label(xC(n),yC(n),1)= i;
        end
    end
    
    %% choose points for convolution one for whole express, one for surface label
    
    %surface label
    num_points_single_cell = num_points(i);
    total_points_Surf = length(cellXs_surface);
    
    %Preallocate to store select points coordinate
    xS = zeros(num_points_single_cell,1);
    yS = xS;
    zS = xS;
    for k=1:num_points_single_cell
        point_number = ceil(rand(1)*total_points_Surf+0.01); %randomly choose a point, ensure point_number not equal to 0
        if point_number>total_points_Surf
            point_number = total_points_Surf;
        end
        try
            xS(k) = ceil(cellXs_surface(point_number)/voxelSize);
            yS(k) = length(yRange)-ceil(cellYs_surface(point_number)/voxelSize);
            zS(k) = ceil(cellZs_surface(point_number)/voxelSize);
            surface_select_points(xS(k),yS(k),zS(k))=surface_select_points(xS(k),yS(k),zS(k))+1;
        catch
           warning('Problem of indexing. 2');
           disp(zS(k));
        end
    end
    
    %whole express
    num_points_single_cell = num_points(i);
    total_points_Whole = length(cellXs);
    
    %Preallocate to store select points coordinate
     xW = zeros(num_points_single_cell,1);
     yW = xW;
     zW = xW;
    for k=1:num_points_single_cell
        point_number = ceil(rand(1)*total_points_Whole+0.01); %randomly choose a point, ensure point_number not equal to 0
        if point_number > total_points_Whole
            point_number = total_points_Whole;
        end
        % make sure xW(k) and so can't be zero for indexing
        try 
            xW(k) = ceil(cellXs(point_number)/voxelSize);
            yW(k) = length(yRange)-ceil(cellYs(point_number)/voxelSize);
            zW(k) = ceil(cellZs(point_number)/voxelSize);
        catch 
            warning('Problem of indexing. 3');
            disp(zW(k));
        end
        try
            whole_select_points(xW(k),yW(k),zW(k))=whole_select_points(xW(k),yW(k),zW(k))+1;
        catch
            warning('Problem of indexing. 4');
        end
    end
end%end for loop
rawLocs=[xLocs,yLocs,zLocs];
end%end function

