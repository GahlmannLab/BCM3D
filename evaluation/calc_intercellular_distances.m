%% calc_intercellular_distances
% work with a single file to test or multiple files
% input is GT or segmentation result
% step 1: get all pixels that are on the boundary of each object
% step 2: iteratively search the minimum distance for each pixel within 
% an object using knnsearch. Get a minimum distance for that object. 
% step 3: Go over all objects to find minimum distances distribution for all objects
% output: meanD_results, medianD_results, and minDistTable for all distances
clear
%% get images (or an image for testing)
addpath(genpath('functions'), genpath('matlab_deskew_1'));
[fileName,dirName] = uigetfile('*','Load the test microscope images',...
    'MultiSelect', 'on');
dataFile = fullfile(dirName, fileName);

%% check and compute
if isequal(dataFile,0)
    error('User cancelled the program');
elseif ischar(dataFile) % only one GT selected
    img = loadtiff(dataFile);
    % get the boundary by first erode and 
    % substract that image from the original image
    surfaceVoxels = img - imerode(img, true(3,3,3)); %surfaceVoxels = img - imerode(img, true(3)); 
    [minDistTable,meanD,medianD] = calc_minIntercellularDist(surfaceVoxels);
else % multiple GT selected
    for i = 1:length(dataFile)
        img = loadtiff(dataFile{i});
        % get the boundary by first erode and 
        % substract that image from the original image
        surfaceVoxels = img - imerode(img, true(3)); 
        [minDistTable,meanD,medianD] = calc_minIntercellularDist(surfaceVoxels);
        meanD_results(i) = meanD;
        medianD_results(i) = medianD;
        minDistTable_results{i} = minDistTable;
    end
    % create and save a table
    finalTable = table(fileName', meanD_results', medianD_results',minDistTable_results', ...
        'VariableNames',{'file names','the mean of minimum intercellular distances',...
        'the median of minimum intercellular distances', 'all_minimum_distances'});
    try
        save('table\finalTable_intercellular_distances_exp','finalTable');
    catch
        mkdir table
        save('table\finalTable_intercellular_distances_exp','finalTable');
    end
end

%% get boundary
function [minDistTable,meanD,medianD] = calc_minIntercellularDist(surfaceVoxels)
    uniqueLabel = unique(surfaceVoxels);
    uniqueLabel = uniqueLabel(uniqueLabel ~= 0);
    gt_index = cell(1,length(uniqueLabel));
    for i=1:length(uniqueLabel) % no zeros
        % calculate density of each cell and add to get total cell volume
        %[x,y,z]= ind2sub(size(img),find(img == uniqueLabel(i)));
        [x,y,z]= ind2sub(size(surfaceVoxels),find(surfaceVoxels == uniqueLabel(i)));
        gt_index{i} =[x y z] ; 
    end

    %% use knnsearch to search for each cell
    % it is based on K-d tree
    % search time is O(log n)
    minDistTable = zeros(1,length(gt_index));
    parfor j = 1:length(gt_index)
        temp_idx = gt_index;
        temp_idx{1,j} = [];
        X = cell2mat(temp_idx');
        Y = gt_index{j};
        [Idx,D] = knnsearch(X,Y);
        minDist = min(D);
        minDistTable(1,j) =  minDist;
        j
    end

    meanD = mean(minDistTable);
    medianD = median(minDistTable);
    fprintf('The median of distance to its nearest neighbor is %f \n', medianD);
    fprintf('The mean of distance to its nearest neighbor is %f \n', meanD );
end
