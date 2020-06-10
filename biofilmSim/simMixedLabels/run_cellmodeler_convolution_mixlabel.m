%Main program to get simulated 3D fluorscence stack of mixlabel biofilm 
%at different mixing ratio, strongly suggest use parallel loop otherwise 
%it will take very long time.
%last updated 06/04/2020
addpath('Function','InputExample');
addpath(genpath('../../utilFiles'));
% Load cellmodeler result and PSF
load('cell_parameter_eachframe.mat');
load('PSF488.mat');
cell_data = parameter_eachcell;
surf_percentage=[0.1,0.3,0.5,0.7,0.9];
% Percentage of membrane labeled cells in the mixlabel biofilms, change
% accordingly, this example will give 5 differnt mixratio mixlabel biofilms
% data.
SBR=8000;% Setting signal to background ratio (SBR)
% Here SBR defined as the total signal of a cell divide the background
% signal
voxelSize = 100; % Setting voxel size
% change outputfolder accordingly
datafolder1 = strcat('Result/mixlabel_raw/');
datafolder2 = strcat('Result/mixlabel_deconv/');
datafolder3 = strcat('Result/surf_gt/');
datafolder4 = strcat('Result/surfinter_gt/');
mkdir(datafolder1);
mkdir(datafolder2);
mkdir(datafolder3);
mkdir(datafolder4);
parfor k=1:length(surf_percentage)
    parforCellmodeler_mixlabel(cell_data, surf_percentage, k, SBR, voxelSize,...
        psf_conv, psf_decon, background, datafolder1, datafolder2, datafolder3, datafolder4);
end