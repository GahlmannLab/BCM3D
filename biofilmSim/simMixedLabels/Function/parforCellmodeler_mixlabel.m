%Function for paralle loop
function parforCellmodeler_mixlabel(cell_data, surf_percentage, k, SBR, voxelSize,...
        psf_conv, psf_decon, background, datafolder1, datafolder2, datafolder3, datafolder4)

s_percent = surf_percentage(k);
percent = num2str(s_percent);
for i=1:length(cell_data)
    index = num2str(i);
    prefix = strcat(index,'_surf_percent_',percent,'_');
    CellModeller_data = cell_data{i};
    Cellmodeler_convolution_mixlabel
end
end