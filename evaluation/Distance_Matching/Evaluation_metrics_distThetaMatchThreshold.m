%% Compute evaluation metrics based on threshold of distance and angle difference
% of the cell major axis
% plot evaluation metrics: Jaccard index, 
% 12-05-2019 Yibo
% can't be cells to be matched twice.
% compute Jaccard counting accuracy instead of DICE
% the inputs are .tif ground truth images and .nii prediction masks.
clear; close all;
addpath(genpath('Distance_Matching'));addpath(genpath('Plot_Metrics_Density'));
path_Label = uigetdir(pwd, 'Select directory of GT labels'); % the directory where saves GT labels
addpath(genpath(path_Label));
%% user input file directory
path_Mask = uigetdir(pwd,'Select directory of 5D inferecned results'); % the directory where saves .nii 5D inferecned results
addpath(genpath(path_Mask));
SLabel = dir(fullfile(path_Label,'*_Label.tif')); % pattern to match ground truth label filenames.
% unzip the files first
SMask = dir(fullfile(path_Mask,'*_out.nii')); % pattern to match output filenames.
if numel(SLabel) ~= numel(SMask)
    error('files numbers dont match.');
end
prompt = {'Output directory'};
dlgtitle = 'Input';
dims = [1 35];
definput = { 'seg_output/folder/'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
output_dir = answer{1};
%% preallocate parameters
density = zeros(length(SLabel), 1);
Density_name = zeros(length(SLabel), 1);
stats = {length(SLabel), 1};
jaccard_counting_accuracy = zeros(length(SLabel), 1);
for k = 1:numel(SLabel)
    %% start processing 
    FLabel = fullfile(path_Label,SLabel(k).name);
    GT_img = open_3D_tiff(FLabel);
    FMask = fullfile(path_Mask,SMask(k).name);
    % second parameter determines the cutoff of the confidance map
    BW=process5Dimages(FMask, 0.94);
    % second parameter determines the degree of img dilation, optional
    instance_seg_img = instance_segmentation_function(BW,2);
    %% somehow the images needed to be rotated to match the correct orientation;
    instance_seg_img = imrotate(instance_seg_img,270);
    instance_seg_img = flip(instance_seg_img,2);
    %% Extract segmented cell information
    seg_cell_stats = regionprops3(instance_seg_img,'all');
    stats{k,1} = seg_cell_stats;
    %% generate a nii file with one integer indicating one instance for testing purpose
    % comment out when not used
    segName = strcat(FLabel,'_seg.nii');
    niftiwrite(instance_seg_img, segName);
    
    %% compute the jaccard Index counting accuracy
    % The input is GT_image, segmentation_image, angle threshold to match (degrees),
    % distance threshold to match (pixels)
    jaccard_counting_accuracy(k,1) = shortest_distance_matching_v3_twoThresh(GT_img, instance_seg_img, 30, 20);
    
    %%compute the density for each GT image
    density(k,1) =  calculate_density(GT_img);
   
end


