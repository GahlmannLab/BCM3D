function [numCells,rod_num, sphere_num, rod_select_points,sphere_select_points,ground_rod_label, ground_sphere_label] = cellVolume_cellmodeler_mixshape(CellModeller_data,voxelSize,simbox,r_percent)
% Arrange each cells within biofilms according to simulated result of
% CellModeller

% Input
% CellModeller_data - simulated results of CellModeller
% simbox - volume size of 3D stack
% voxelSize - voxel size
% r_percent - percentage of rod shaped cells

% Output
% numCells - number of cells
% rod_num - number of rod shaped cells
% sphere_num - number of sphere shaped cells
% rod_select_points - selected fluorophores in rod shaped cells
% sphere_select_points - selected fluorophores in sphere shaped cells
% ground_rod_label - ground truth of rod shaped cells
% ground_sphere_label - ground truth of sphere shaped cells

numCells = size(CellModeller_data,1);
sp = 30; % grid unit size to arrange individual cells
% translation matrix
COM = CellModeller_data(:, 4:6) + 0.5 * [simbox(1),simbox(2),1000];
% cell orientation
cell_direction_vectors = CellModeller_data(:, 7:9);

%Creating voxel grid
xRange = 0:voxelSize:simbox(1);
yRange = 0:voxelSize:simbox(2);
zRange = 0:voxelSize:simbox(3);

%preallocate
ground_rod_label = zeros(length(xRange),length(yRange),length(zRange));
ground_sphere_label = ground_rod_label;
rod_select_points = zeros(length(xRange),length(yRange),length(zRange));
sphere_select_points = rod_select_points;

%Determine whether a cell is rod shape or sphere shape
rod_num = round(numCells*r_percent);
sphere_num = numCells-rod_num;
index = randperm(numCells);
rod_index = index(1:rod_num);
sphere_index = index(rod_num+1:end);

%number of fluorophores chosen from each cell for convolution
num_points=zeros(numCells, 1);
for i=1:numCells
    num_points(i)=round(500+rand(1)*500);
end

%% rod shaped cells
if rod_num > 0
lc = CellModeller_data(:, 3);% cell length
for i = 1:rod_num
    j = rod_index(i);
    dx = COM(j,1);
    dy = COM(j,2);
    dz = COM(j,3);
    
    %% generate points for rod shpaed cells
    r = 450;% radius of rod shpaed cells
    [All_points,~] = CellLocsCaps_CellModeler(r,lc(j),sp);
    
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
    
    
    %% Converting continuous cell locs to voxel coordinates and filling ground
    xC = zeros(length(cellXs),1);
    yC = xC;
    zC = xC;
    for n=1:length(cellXs)
        xC(n) = ceil(cellXs(n)/voxelSize);
        yC(n) = length(yRange)-ceil(cellYs(n)/voxelSize); % 
        zC(n) = ceil(cellZs(n)/voxelSize);
        ground_rod_label(xC(n),yC(n),zC(n))= i;
    end
    
    %% Choose fluorophores
    
    num_points_single_cell = num_points(j);
    total_points = length(cellXs);
    
    %Preallocate to store selected fluorophores coordinate
    xS = zeros(num_points_single_cell,1);
    yS = xS;
    zS = xS;
    %randomly choose fluorophores
    for k=1:num_points_single_cell
        point_number = ceil(rand(1)*total_points+0.01); 
        if point_number>total_points
            point_number = total_points;
        end
        xS(k) = round(cellXs(point_number)/voxelSize);
        yS(k) = length(yRange)-round(cellYs(point_number)/voxelSize);
        zS(k) = round(cellZs(point_number)/voxelSize);
        rod_select_points(xS(k),yS(k),zS(k))=rod_select_points(xS(k),yS(k),zS(k))+1;
    end
end % end rod shpae loop
end % end if 

%% Sphere shape cell
if sphere_num > 0
for i = 1:sphere_num
    j = sphere_index(i);
    dx = COM(j,1);
    dy = COM(j,2);
    dz = COM(j,3);
    
    %% generate points for the cell
    % need transfer size of rod shape to sphere shape cells, here CellModeller
    % only simulate rod shpaed cells
    lspher = rand(1)*100+50;
    r = 450;%radius of sphere shaped cells
    [All_points,~] = CellLocsCaps_CellModeler(r,lspher,sp);
    
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
    

    %% Converting continuous cell locs to voxel coordinates and filling ground
    xC = zeros(length(cellXs),1);
    yC = xC;
    zC = xC;
    for n=1:length(cellXs)
        xC(n) = ceil(cellXs(n)/voxelSize);
        yC(n) = length(yRange)-ceil(cellYs(n)/voxelSize); % 
        zC(n) = ceil(cellZs(n)/voxelSize);
        ground_sphere_label(xC(n),yC(n),zC(n))= i;
    end
   %% Choose fluorophores for convolution

    total_points = length(cellXs);
    
    %Preallocate to store select fluorophores coordinate
    xS = zeros(num_points_single_cell,1);
    yS = xS;
    zS = xS;
    %randomly choose fluorophores
    for k=1:num_points_single_cell
        point_number = ceil(rand(1)*total_points+0.01); 
        if point_number>total_points
            point_number = total_points;
        end
        xS(k) = round(cellXs(point_number)/voxelSize);
        yS(k) = length(yRange)-round(cellYs(point_number)/voxelSize);
        zS(k) = round(cellZs(point_number)/voxelSize);
        sphere_select_points(xS(k),yS(k),zS(k))=sphere_select_points(xS(k),yS(k),zS(k))+1;
    end

end%end sphere for loop
end% end if 


end%end function
