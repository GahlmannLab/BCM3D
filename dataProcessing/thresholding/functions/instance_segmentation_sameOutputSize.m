
%% Instance segmentation for 3D Lattice light sheet microscopy
% Load the inferenced data from deep learning model(niftinet/unet)
% Usually input is .nii
% the output is instance segmented image saved in .nii 
% 
%%%

%%Load the file and separate each objects
function [segmented_img] = instance_segmentation_sameOutputSize(BW, filterSize,varargin)

% [BW]=process5Dimages(fileName);
%% use bwareaopen to filter big and small cells
% bwconncomp to get connected cells
seg_img = bwareaopen(BW,filterSize,6);
CC = bwconncomp(seg_img,6);
stats = regionprops3(CC,'VoxelList');
cells = stats.VoxelList;
% save the parameter, size of the image
sz = size (seg_img);
new_cell = cells;
%% Instance segmentation in matlab
m=zeros(sz);
m = uint16(m);
% Transform coordinates into the matrix where each cell is labeled as an
% integer number
for cell_counter = 1:length(new_cell)
    each_cell = new_cell(cell_counter,:);
    each_cell = cell2mat(each_cell);
    each_cell = int16(each_cell);
    
    Y = each_cell(:,1);
    X = each_cell(:,2);
    Z = each_cell(:,3);
    %Transform coordinates into the matrix where each 
    %cell is labeled as an integer as cell_counter
    for i = 1:length(X)
        m((X(i)),Y(i),Z(i)) = cell_counter;
    end
end
%% Dilation of each cell since cell boudaries have been eroded controlled by varargin
if  nargin == 2
    segmented_img = m;
elseif  nargin == 3 
    if varargin{1} ~= 0
        SE = strel('sphere',varargin{1}); 
        segmented_img = imdilate(m,SE,'same');
        segmented_img = segmented_img(1:sz(1),1:sz(2),1:sz(3));
    else
       segmented_img = m;
    end
else
    error('second parameter has to be a positive interger');
end


