%Function for paralle loop
function parforCellmodeler_mixlabel(cell_data, surf_percentage, k)
s_percent = surf_percentage(k);
percent = num2str(s_percent);
% Change output folder name accordingly
datafolder1 = strcat('Result/mixlabel_raw/');
datafolder2 = strcat('Result/mixlabel_deconv/');
datafolder3 = strcat('Result/surf_gt/');
datafolder4 = strcat('Result/surfinter_gt/');
mkdir(datafolder1);
mkdir(datafolder2);
mkdir(datafolder3);
mkdir(datafolder4);
for i=1:length(cell_data)
    index = num2str(i);
    prefix = strcat(index,'_surf_percent_',percent,'_');
    CellModeller_data = cell_data{i};
    SBR=8000;% Setting signal to background ratio (SBR)
    voxelSize = 100; % Voxel size
    Cellmodeler_convolution_mixlabel
end
end