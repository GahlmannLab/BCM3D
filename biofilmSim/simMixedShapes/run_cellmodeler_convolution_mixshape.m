%Main program to get simulated 3D fluorscence stack of mix shape biofilms 
%at different mixing ratio, strongly suggest use parallel loop otherwise it 
%will take very long time.
%last updated 05/22/2020
addpath('Function','InputExample');
addpath(genpath('../../utilFiles'));
load('cell_parameter_eachframe.mat');
cell_data = parameter_eachcell;
rod_percentage= [0.1,0.3,0.5,0.7,0.9];
% Percentage of rod shaped cells in the mixshape biofilms, change
% accordingly, this example will give 5 differnt mixratio mixshape biofilms
% data.
parfor k=1:length(rod_percentage)
    parforCellmodeler_mixshape(cell_data, rod_percentage, k);
end