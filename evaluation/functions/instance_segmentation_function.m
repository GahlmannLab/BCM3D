
%% Instance segmentation for 3D Lattice light sheet microscopy
% Load the inferenced data from deep learning model(niftinet/unet)
% Usually input is .nii
% the output is instance segmented image saved in .nii 
% 
% updated: 20200113 YW
% compatible with niftiread_v2/niftiwrite_v2
%%%

%%Load the file and separate each objects
function [segmented_img] = instance_segmentation_function(BW, varargin)

% [BW]=process5Dimages(fileName);
%% use bwareaopen to filter big and small cells
% bwconncomp to get connected cells
seg_img = bwareaopen(BW,50,6);
CC = bwconncomp(seg_img,6);
stats = regionprops3(CC,'VoxelList');
% stats = regionprops3(CC,'VoxelIdxList');
cells = stats.VoxelList;
% cells = stats.VoxelIdxList;
% save the parameter, size of the image
sz = size (seg_img);
%% select and plot on the undersegmented cells that we want to check
% BW = zeros(sz);
% for cell = 1:length(cells)
%     if stats.ConvexVolume(cell,1) > 2000 && stats.SurfaceArea(cell,1) > 1600
%         a = int16(cell2mat(cells(cell,1)));
%         X = a(:,1);
%         Y = a(:,2);
%         Z = a(:,3);
%         %plot3(X,Y,Z,'.');
%         %hold on
%         for i = 1:length(X)
%             BW((X(i)),Y(i),Z(i)) = 1;
%         end
%         cells{cell,1} = [];
%     end
% end
% %%
% % erode cell and add back, sometimes can decrease over-segmentation
% SE = strel('sphere',2);
% BW1 = imerode(BW, SE);
% BW1_open = bwareaopen(BW1,60,6);
% % filter agian out big and small cells
% BW1_open = xor(BW1_open , bwareaopen(BW1_open , 5000)); 
% 
% % generate new_cell which is a matlab cell array for future use
% CC_cluster = bwconncomp(BW1_open,6);
% stats_cluster = regionprops3(CC_cluster,'all');
% cells_cluster = stats_cluster.VoxelList;
% new_cell = cat(1,cells,cells_cluster);
% new_cell = new_cell(~cellfun('isempty',new_cell));
%%
% plot new_cell as new figure, to see if it can match with
%Data_visualization file 
% figure;
% hold on;axis equal;
% cmap = jet(size(new_cell,2));
% for cell = 1:length(new_cell);
%         colorm = rand(1,3);
%         a = int16(cell2mat(new_cell(cell,1)));
%         X = a(:,1);
%         Y = a(:,2);
%         Z = a(:,3);
%         plot3(X,Y,Z,'.','color',colorm,'MarkerSize',10);
%         axis([0 300 0 300 ]);
%         hold on
% end
%% Use plot3 function to plot 3D graphs of instance segmentation result

%%
% figure;
% for cell = 1:length(new_cell)
%     
%     a = cell2mat(new_cell(cell,1));
%     X = a(:,1);
%     Y = a(:,2);
%     Z = a(:,3);
%     plot3(X,Y,Z,'.');
%     hold on
% end
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
    
    X = each_cell(:,1);
    Y = each_cell(:,2);
    Z = each_cell(:,3);
    %Transform coordinates into the matrix where each 
    %cell is labeled as an integer as cell_counter
    for i = 1:length(X)
        m((Y(i)),X(i),Z(i)) = cell_counter; % this is compatible with tiff and nifti_v2
    end
%     m(cells{cell_counter}) = cell_counter;

end
%% Dilation of each cell since cell boudaries have been eroded controlled by varargin
if  nargin == 1
    segmented_img = m;
elseif  nargin == 2
    SE = strel('sphere',varargin{1}); 
    segmented_img = imdilate(m,SE);
else
    error('second parameter has to be a positive interger');
end


