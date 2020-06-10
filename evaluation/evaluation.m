%% Compute evaluation metrics
% plot evaluation metrics: Jaccard index, Dice, False positive percentage,
% False negative percentage, Distance of major axis, and relative angle of
% major axis. Counting accuracy added (20190913)
%  after each match the GT cell gets deleted. 
% There can't be cells that are matched twice.
% compute Jaccard counting accuracy instead of DICE
% the inputs are .tif ground truth images and segmentation masks. 
% change nifti to nifti_v2 to avoid flipping and rotation
% updated 20200113 YW
clear; close all;
addpath(genpath('functions'));addpath(genpath('../utilFiles'));
path_Label = uigetdir(pwd, 'Select directory of GT labels'); % the directory where saves GT labels
addpath(genpath(path_Label));
%% user input file directory
path_Mask = uigetdir(pwd,'Select directory of segmentation results'); % the directory where saves .nii 5D inferecned results
addpath(genpath(path_Mask));
SLabel = dir(fullfile(path_Label,'*.tif')); % pattern to match ground truth label filenames.
% unzip the files first
SMask = dir(fullfile(path_Mask,'*.tif')); % pattern to match output filenames.
if numel(SLabel) ~= numel(SMask)
    error('files numbers dont match.');
end
%% user input name of the saved file and directory 
prompt = {'Saved name'};
dlgtitle = 'Input';
dims = [1 60];
definput = {'Evaluation_metrics.mat'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
eval_name = answer{1};
%% preallocate parameters
Dice_eachFrame = zeros(length(SLabel), 1);
Jaccard_eachFrame = zeros(length(SLabel), 1);
match_percent = zeros(length(SLabel), 1);
density = zeros(length(SLabel), 1);
FN_eachFrame =  zeros(length(SLabel), 1);
FP_eachFrame =  zeros(length(SLabel), 1);
distance_distribution_eachframe = {length(SLabel), 1};
angle_distribution_eachframe = {length(SLabel), 1};
figure_names = zeros(length(SLabel), 1);
SBR_name = zeros(length(SLabel), 1);
Density_name = zeros(length(SLabel), 1);
stats = {length(SLabel), 1};
jaccard_counting_accuracy = zeros(length(SLabel), 1);
dice_counting_accuracy = zeros(length(SLabel), 1);
parfor k = 1:numel(SLabel)
    FLabel = fullfile(path_Label,SLabel(k).name);
    seg_img = fullfile(path_Mask,SMask(k).name); 
    figure_name = replace(SLabel(k).name, '_deconvolvedwholeexp_Label.tif', '');
    %% find the string name of SBR and density for plotting
    index = strfind(figure_name, '_');
    SBR_name(k,1) = str2double(figure_name(index(2)+1: index(3)-1));
    Density_name(k,1) =str2double( figure_name(index(4)+1: index(5)-1));
    figure_names(k,1) = str2double(figure_name(index(1)-1));
    %% start processing 
    GT_img = loadtiff(FLabel);
    % second parameter determines the degree of img dilation, optional
    instance_seg_img = loadtiff(seg_img);
    %% Extract segmented cell information
    seg_cell_stats = regionprops3(instance_seg_img,'all');
    stats{k,1} = seg_cell_stats;
    % compute the jaccard Index,
    % Dice, False positive percentage,
    % False negative percentage,
    [Jaccard, Dice, match_num, jaccardAccuracy, diceAccuracy, FN, FP] = jaccard_Index_3D_Evaluation_v4_deleteMatched_JaccardAccuracy(GT_img, instance_seg_img, 0.5)
    Jaccard_eachFrame(k,1)= Jaccard;
    Dice_eachFrame(k,1)= Dice;
    FN_eachFrame(k,1)= FN;
    FP_eachFrame(k,1)= FP;
    match_percent(k,1) = match_num;
    jaccard_counting_accuracy(k,1) = jaccardAccuracy;
    dice_counting_accuracy(k,1) = diceAccuracy;
    %%
    %compute the density for each GT image
    tile = [64 64 8]; %tile size
    incrementxy = 8; %xy overlaps
    incrementz = 4; % z overlaps
    [~, density(k,1)] =  calc_tile_maximum_density(GT_img,tile,incrementxy,incrementz);
    %%
    % Plot the Jaccard index vs. cell density   
    message1 = ['the Jaccard Index is ', num2str(Jaccard)];
    message2 = ['the Jaccard counting accuracy is ', num2str(jaccardAccuracy)];
    message3 = ['the Dice counting accuracy is ', num2str(diceAccuracy)];
    disp(message1);
    disp(message2);
    disp(message3);
end
distance_distribution_eachframe=distance_distribution_eachframe(:,1);
angle_distribution_eachframe=angle_distribution_eachframe(:,1);
%% create a final table to summarize all the metrics
Evaluation_metrics = table(figure_names,SBR_name, Density_name,...
    Dice_eachFrame, Jaccard_eachFrame, match_percent,...
    density, FN_eachFrame, jaccard_counting_accuracy, FP_eachFrame);
% save the table as .mat in the matFileOutput folder
try
    full_eval_name = strcat('matFileOutput/',eval_name );
    save(full_eval_name, 'Evaluation_metrics');
catch 
    mkdir matFileOutput
    full_eval_name = strcat('matFileOutput/',eval_name );
    save(full_eval_name, 'Evaluation_metrics');
end

