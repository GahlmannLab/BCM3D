#### Simulate 3D fluorescence mixlabeled biofilms stacks and their corresponding GTs.**

###### Mixlabeled biofilms means multi population biofilms that cells have different labelling protocols ( e.g. membrane  vs. cytosolic&membrane labeling).



The master script is run_cellmodeler_convolution_mixlabel.m .

*inputs:*

It reads CellModeller simulation results from cell_parameter_eachframe.mat.

It uses parforCellmodeler_mixlabel function for parallel computing, change parameter accordingly in parforCellmodeler_mixlabel.

*outputs:*
It saves .tif ground truth separately into Result/surf_gt and Result/surfinter_gt

It saves raw simulated fluorescence stack .tif into Result/mixlabel_raw

It saves deconvolved simulated fluorescence stack .nii into Result/mixlabel_deconv

