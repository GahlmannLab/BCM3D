function [numCells,surf_num, surf_interior_num,surf_select_points,surf_interior_select_points,ground_surf_label, ground_surf_interior_label] = cellVolume_cellmodeler_mixlabel(CellModeller_data,voxelSize,simbox,s_percent)
% Arrange each cells within biofilms according to simulated result of
% CellModeller

% Input
% CellModeller_data - simulated results of CellModeller
% simbox - volume size of 3D stack
% voxelSize - voxel size
% r_percent - percentage of rod shaped cells

% Output
% numCells - number of cells
% surf_num - number of membrane labeled cells
% surf_interior_num - number of membrane&cytosolic labeled cells
% surf_select_points - selected fluorophores in membrane labeled cells
% surf_interior_select_points - selected fluorophores in membrane&cytosolic labeled cells
% ground_surf_label - ground truth of membrane labeled cells
% ground_surf_interior_label - ground truth of membrane&cytosolic labeled cells


numCells = size(CellModeller_data,1);
sp = 30; % grid unit size to arrange individual cells

% translation matrix
COM = CellModeller_data(:, 4:6) + 0.5 * [simbox(1),simbox(2),1000];
% cell orientation
cell_direction_vectors = CellModeller_data(:, 7:9);
% cell radius and length
r = CellModeller_data(:, 2);
lc = CellModeller_data(:, 3);


%% Creating voxel grid 
xRange = 0:voxelSize:simbox(1);
yRange = 0:voxelSize:simbox(2);
zRange = 0:voxelSize:simbox(3);

%preallocate
ground = zeros(length(xRange),length(yRange),length(zRange));
ground_surf_label = ground;
ground_surf_interior_label = ground;
surf_select_points = zeros(length(xRange),length(yRange),length(zRange));
surface_select_points_temp = surf_select_points;
interior_select_points_temp = surf_select_points;

%Determine whether a cell is membrane labeled or membrane&cytosolic labeled
det = rand(1,numCells);
surf_num = sum(det<s_percent);
surf_interior_num = sum(det>=s_percent);
surf_index = find(det<s_percent);
surf_interior_index = find(det>=s_percent);

%Number of fluorophores chosen from each cell
num_points=zeros(numCells, 1);
for i=1:numCells
    num_points(i)=round(500+rand(1)*500);
end

%% Membrane labeled cells
if surf_num > 0
for i = 1:surf_num
    j = surf_index(i);
    dx = COM(j,1);
    dy = COM(j,2);
    dz = COM(j,3);
    
    %% generate points for the cell
    [All_points,surface_points] = CellLocsCaps_CellModeler(r(j),lc(j),sp);
    
    original_unit_vector = [1 0 0];
    final_unit_vector = cell_direction_vectors(j, :);
    rotated_cell_points = rotCellWithVector(All_points, original_unit_vector, final_unit_vector);
    
    cellXs = rotated_cell_points(:,1);
    cellYs = rotated_cell_points(:,2);
    cellZs = rotated_cell_points(:,3);
    
    % move the cell
    cellXs = cellXs + dx;
    cellYs = cellYs + dy;
    cellZs = cellZs + dz;
    
    
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
    % Preallocate
    xC = zeros(length(cellXs),1);
    yC = xC;
    zC = xC;
    for n=1:length(cellXs)
        xC(n) = ceil(cellXs(n)/voxelSize);
        yC(n) = length(yRange)-ceil(cellYs(n)/voxelSize); % 
        zC(n) = ceil(cellZs(n)/voxelSize);
        ground(xC(n),yC(n),zC(n))= 1;
        ground_surf_label(xC(n),yC(n),zC(n))= i;
    end
    
    %% choose fluorophores
    
    num_points_single_cell = num_points(j);
    total_points_Surf = length(cellXs_surface);
    
    %Preallocate to store select points coordinate
    xS = zeros(num_points_single_cell,1);
    yS = xS;
    zS = xS;
   %randomly choose fluorphores
    for k=1:num_points_single_cell
        point_number = ceil(rand(1)*total_points_Surf+0.01);
        if point_number>total_points_Surf
            point_number = total_points_Surf;
        end
        xS(k) = ceil(cellXs_surface(point_number)/voxelSize);
        yS(k) = length(yRange)-ceil(cellYs_surface(point_number)/voxelSize);
        zS(k) = ceil(cellZs_surface(point_number)/voxelSize);
        surf_select_points(xS(k),yS(k),zS(k))=surf_select_points(xS(k),yS(k),zS(k))+1;
    end
end % end membrane labeled cells loop
end % end if 

%% Membrane&cytosolic labeled cells
if surf_interior_num > 0
for i = 1:surf_interior_num
    j = surf_interior_index(i);
    dx = COM(j,1);
    dy = COM(j,2);
    dz = COM(j,3);
    
    %% generate points for the cell
    [All_points,surface_points] = CellLocsCaps_CellModeler(r(j),lc(j),sp);
    
    original_unit_vector = [1 0 0];
    final_unit_vector = cell_direction_vectors(j, :);
    rotated_cell_points = rotCellWithVector(All_points, original_unit_vector, final_unit_vector);
    
    cellXs = rotated_cell_points(:,1);
    cellYs = rotated_cell_points(:,2);
    cellZs = rotated_cell_points(:,3);
    
    % move the cell
    cellXs = cellXs + dx;
    cellYs = cellYs + dy;
    cellZs = cellZs + dz;
    
    %% generate membrane points for the cell
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
        ground(xC(n),yC(n),zC(n))= 1;
        ground_surf_interior_label(xC(n),yC(n),zC(n))= i;
    end
    
    %% choose fluorophores
    
    num_points_single_cell = num_points(j);
    total_points_Surf = length(cellXs_surface);
    %membrane
    %Preallocate to store selected fluorophores coordinate
    xS = zeros(num_points_single_cell,1);
    yS = xS;
    zS = xS;
    %randomly choose fluorophores
    for k=1:num_points_single_cell
        point_number = ceil(rand(1)*total_points_Surf+0.01); 
        if point_number>total_points_Surf
            point_number = total_points_Surf;
        end
        xS(k) = ceil(cellXs_surface(point_number)/voxelSize);
        yS(k) = length(yRange)-ceil(cellYs_surface(point_number)/voxelSize);
        zS(k) = ceil(cellZs_surface(point_number)/voxelSize);
        surface_select_points_temp(xS(k),yS(k),zS(k))=surface_select_points_temp(xS(k),yS(k),zS(k))+1;
    end
    %cytosolic
    %Preallocate to store selected fluorophores coordinate
    total_points_Interior = length(cellXs);
     xW = zeros(num_points_single_cell,1);
     yW = xW;
     zW = xW;
    for k=1:num_points_single_cell
        point_number = ceil(rand(1)*total_points_Interior+0.01); 
        if point_number > total_points_Interior
            point_number = total_points_Interior;
        end
        xW(k) = ceil(cellXs(point_number)/voxelSize);
        yW(k) = length(yRange)-ceil(cellYs(point_number)/voxelSize);
        zW(k) = ceil(cellZs(point_number)/voxelSize);
        interior_select_points_temp(xW(k),yW(k),zW(k))=interior_select_points_temp(xW(k),yW(k),zW(k))+1;
    end
    %sum membrane and cytosolic fluorophores
    surf_interior_select_points = interior_select_points_temp + surface_select_points_temp;
end%end membrane&cytosolic labeled cells for loop
else
    surf_interior_select_points = zeros(length(xRange),length(yRange),length(zRange));
end% end if 

end%end function
