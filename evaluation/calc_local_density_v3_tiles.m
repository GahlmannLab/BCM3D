%% test_local_density_v3_multipleFiles
% work with a single file to test or multiple files
% tile GT images and calculate local density
% local density = cell volume / tile volume
% two output metrics: maximum local density and the mean of top 10 local densities
clear
%% tile parameters
%change as needed
saveName = 'mixshape_density';
tile = [64 64 8]; %tile size
incrementxy = 8; %xy overlaps
incrementz = 4; % z overlaps

%% select data
addpath(genpath('functions'), genpath('matlab_deskew_1'));
[fileName,dirName] = uigetfile('*','Load the test microscope images', ...
    'MultiSelect', 'on');
dataFile = fullfile(dirName, fileName);

%% check and compute
if isequal(dataFile,0)
    error('User cancelled the program');
elseif ischar(dataFile) % only one GT selected
    img = loadtiff(dataFile);
    [maxLocalDensity, mean10LocalDensity] = calc_tile_maximum_density(img,tile,incrementxy,incrementz);
else % multiple files are selected, save as a table
    % calc density of bilfilms
    mean10LocalDensity_results = [];
    maxLocalDensity_results = [];
    parfor i = 1:length(dataFile)
        img = loadtiff(dataFile{i});
        [maxLocalDensity, mean10LocalDensity] = calc_tile_maximum_density(img,tile,incrementxy,incrementz);
        mean10LocalDensity_results(i) = mean10LocalDensity;
        maxLocalDensity_results(i) = maxLocalDensity;
    end
    % create and save a table
    finalTable = table(fileName', mean10LocalDensity_results', maxLocalDensity_results',...
        'VariableNames',{'file names','mean10LocalDensity','maxLocalDensity'});
    try
        save(['table\',saveName],'finalTable');
    catch
        mkdir table
        save(['table\',saveName],'finalTable');
    end
end

%% devide image into local cubes 
function [maxLocalDensity, mean10LocalDensity] = calc_tile_maximum_density(img,tile,incrementxy,incrementz )
%% devide image into local cubes 
    %parameters
%     tile = [64 64 8];
%     incrementxy = 8; %xy overlaps
%     incrementz = 4; % z overlaps

    kernalVolume = tile(1) * tile(2) * tile(3);
    blockSizeRow = tile(1);blockSizeCol = tile(2);blockSizeZ = tile(3);

    %% cal density of each block
    % blockSizeRow: Rows in block
    % blockSizeCol: Columns in block
    [nrows, ncols, nz] = size(img);

    % Calculate size of each block by rows and columns
    n = 1;
    for i=1:incrementxy:nrows-blockSizeRow+1
         for j = 1:incrementxy:ncols-blockSizeCol+1
             for k = 1:incrementz:nz-blockSizeZ+1
                tileImages = img(i:i+blockSizeRow-1,...
                    j:j+blockSizeCol-1,...
                    k:k+blockSizeZ-1); 

                 if sum(tileImages,'all')== 0
                    continue
                    else   
                        cellVolume = get_cellVolume(tileImages);
                        density(n,1) =  cellVolume/kernalVolume * 100;
                 end
                 n = n +1;
             end
          end
    end

    %density_noOutlier = rmoutliers(density(density~= 0));
    maxLocalDensity = max(density);
    mean5LocalDensity = mean(maxk(density,5));
    mean10LocalDensity = mean(maxk(density,10));
    histogram(density(density ~= 0),20)
    xlabel('local density (%)')
    ylabel('counts')
    fprintf('The maximum of local density of the selected image is %f percent \n', maxLocalDensity);
    fprintf('The mean of the maximum local density of the top five tiles is %f percent \n', mean5LocalDensity);
    fprintf('The mean of the maximum local density of the top ten tiles is %f percent \n', mean10LocalDensity);
end
%% calculate density of the biofilm volume
% this is test_density_v2
% it work with any type of labels
function cellVolume = get_cellVolume(GT_img)
    sz = size(GT_img);
    l_gt = unique(GT_img);
    l_gt(l_gt == 0)=[];
    gt_index = cell(1,length(l_gt));
    cellVolume = 0;
    for i=1:length(l_gt)
        % calculate density of each cell and add to get total cell volume
        r1 = find(GT_img == l_gt(i));
        gt_index{i} = r1; 
        cellVolume = cellVolume+ length(r1);
    end
end



