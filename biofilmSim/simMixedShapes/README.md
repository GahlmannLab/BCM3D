#### Simulate 3D fluorescence mixshaped biofilms stacks and their corresponding GTs.**

###### Mixshaped biofilms means multi population biofilms that cells have different shapes (e.g. spherical vs. rod-shaped cells). 

The master script is run_cellmodeler_convolution_mixshape.m .

*inputs:*

It reads CellModeller simulation results from cell_parameter_eachframe.mat.

It uses parforCellmodeler_mixlabel function for parallel computing, change parameter accordingly in parforCellmodeler_mixlabel.

Here, spherical shaped cells are not simulated by CellModeller directly, cell length for spherical shaped cells are adjusted in cellVolume_cellmodeller.

*outputs:*
It saves .tif ground truth separately into Result/rod_gt and Result/sphere_gt

It saves raw simulated fluorescence stack .tif into Result/mixshape_raw

It saves deconvolved simulated fluorescence stack .tif into Result/mixshape_deconv