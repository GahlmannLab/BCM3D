%Function for paralle loop
function parforCellmodeler_mixshape(cell_data, rod_percentage, k)
r_percent = rod_percentage(k);
percent = num2str(r_percent);
% Change output folder name accordingly
datafolder1 = strcat('Result/mixshape_raw/');
datafolder2 = strcat('Result/mixshape_deconv/');
datafolder3 = strcat('Result/rod_gt/');
datafolder4 = strcat('Result/sphere_gt/');
mkdir(datafolder1);
mkdir(datafolder2);
mkdir(datafolder3);
mkdir(datafolder4);
for i=1:length(cell_data)
    index = num2str(i);
    prefix = strcat(index,'_rod_percent_',percent,'_');
    CellModeller_data = cell_data{i};
    SBR=8000;% Setting signal to background ratio (SBR)
    voxelSize = 100; % Voxel size
    Cellmodeler_convolution_mixshape
end
end