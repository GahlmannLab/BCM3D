%Main program to get simulated 3D fluorscence stack of mixlabel biofilm 
%at different mixing ratio, strongly suggest use parallel loop otherwise 
%it will take very long time.
%last updated 05/22/2020
addpath('Function','InputExample');
addpath(genpath('../../utilFiles'));
load('cell_parameter_eachframe.mat');
cell_data = parameter_eachcell;
surf_percentage=[0.1,0.3,0.5,0.7,0.9];
% Percentage of membrane labeled cells in the mixshape biofilms, change
% accordingly, this example will give 5 differnt mixratio mixshape biofilms
% data.
parfor k=1:length(surf_percentage)
    parforCellmodeler_mixlabel(cell_data, surf_percentage, k);
end