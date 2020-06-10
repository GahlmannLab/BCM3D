%Main program to get simulated 3D fluorscence stack of mix shape biofilms 
%at different mixing ratio, strongly suggest use parallel loop otherwise it 
%will take very long time.
%last updated 05/22/2020
addpath('Function','InputExample');
addpath(genpath('../../utilFiles'));
% Load cellmodeler result and PSF
load('cell_parameter_eachframe.mat');
load('PSF488.mat');
cell_data = parameter_eachcell;
rod_percentage=[0.1,0.3,0.5,0.7,0.9];
% Percentage of rod cells in the mixshape biofilms, change
% accordingly, this example will give 5 differnt mixratio mixshape biofilms
% data.
SBR=8000;% Setting signal to background ratio (SBR)
% Here SBR defined as the total signal of a cell divide the background
% signal
voxelSize = 100; % Setting voxel size
% change outputfolder accordingly
% Change output folder name accordingly
datafolder1 = strcat('Result/mixshape_raw/');
datafolder2 = strcat('Result/mixshape_deconv/');
datafolder3 = strcat('Result/rod_gt/');
datafolder4 = strcat('Result/sphere_gt/');
mkdir(datafolder1);
mkdir(datafolder2);
mkdir(datafolder3);
mkdir(datafolder4);
parfor k=1:length(rod_percentage)
    parforCellmodeler_mixshape(cell_data, rod_percentage, k, SBR, voxelSize,...
        psf_conv, psf_decon, background, datafolder1, datafolder2, datafolder3, datafolder4);
end