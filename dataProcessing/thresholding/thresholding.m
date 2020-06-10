% Instance segmentation of multiple files
% input .gz zipped nii files; output segmented nii files.
% upToDate: 20190923 Yibo
% required function: instance_segmentation_sameOutputSize
clear;
%% user input of files and parameters
addpath(genpath('functions'));
[files, path]= uigetfile('*','multiselect','on');
destDirectory = strcat(path, 'instance_seg_result/');
% Enter postprocessing parameters
prompt = {'Which channel to process:', 'Threshold value',  'number of pixel to do dilation'};
dlgtitle = 'Enter postprocessing parameters';
definput = {'2', '0.95',  '1'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,[1 40],definput,opts);
channel = str2double(answer{1});thresh = str2double(answer{2}); 
dilationSize = str2double(answer{3});
%% make output folder
try 
    mkdir (destDirectory);
catch 
    warning('directories exist');
end
%% check how many input images and process
if  iscell(files) == 1 %% means multiple images are processed
    for i = 1: length(files)
        try
            unzipFile = gunzip(strcat(path,files{1,i}));
        %change the second parameter which indicates class/channel of the output image.
            BW = process5Dimages_multiclass(unzipFile{1}, thresh, channel); 
             % new name tif
            name = erase(files{1,i}, '.nii.gz');
        % second parameter indicates the size minimum of an object, the third parameter indicates dilation.
        catch % no need to unzip
            BW = process5Dimages_multiclass(strcat(path,files{1,i}), thresh, channel);
            name = erase(files{1,i}, '.nii');
        end
        
        segmented_img = instance_segmentation_function(BW, dilationSize);
        % somehow the images needed to be rotated to match the correct orientation;
        newName =  strcat(destDirectory, name, '_thresh', answer{2},...
            '_dilate',answer{3},  '_segmented.tif');
        % write tif
        write3Dtiff_V2(segmented_img,newName);
    end
else  % process only one image
    try
        unzipFile = gunzip(strcat(path,files));
        BW = process5Dimages_multiclass(unzipFile{1}, thresh, channel); 
        name = erase(files, '.nii.gz');
    catch
         BW = process5Dimages_multiclass(strcat(path,files), thresh, channel); 
        name = erase(files, '.nii');
    end
    segmented_img = instance_segmentation_function(BW, dilationSize);
    newName =  strcat(destDirectory, name, '_thresh', answer{2},...
            '_dilate',answer{3},  '_segmented.tif');
    % write tif
    write3Dtiff_V2(segmented_img,newName);
end

%% optional plot
%plotIsoSurface(segmented_img);
