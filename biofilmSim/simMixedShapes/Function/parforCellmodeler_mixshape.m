%Function for paralle loop
function parforCellmodeler_mixshape(cell_data, rod_percentage, k, SBR, voxelSize,...
        psf_conv, psf_decon, background, datafolder1, datafolder2, datafolder3, datafolder4)

r_percent = rod_percentage(k);
percent = num2str(r_percent);
for i=1:length(cell_data)
    index = num2str(i);
    prefix = strcat(index,'_rod_percent_',percent,'_');
    CellModeller_data = cell_data{i};
    Cellmodeler_convolution_mixshape
end
end