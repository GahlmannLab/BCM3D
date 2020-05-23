%% Compute evaluation metrics
% plot evaluation metrics: Jaccard index, Dice, False positive percentage,
% False negative percentage, Distance of major axis, and relative angle of
% major axis. Counting accuracy added (20190913)
%  after each match the GT cell gets deleted. There
% can't be cells to be matched twice.
% compute Jaccard counting accuracy instead of DICE
% the inputs are .tif ground truth images and .nii prediction masks. Updated 11222019 YW
% change nifti to nifti_v2 to avoid flipping and rotation
% updated 20200113 YW
clear; close all;
addpath(genpath('Distance_Matching'));addpath(genpath('functions'));addpath(genpath('../utilFiles'));
path_Label = uigetdir(pwd, 'Select directory of GT labels'); % the directory where saves GT labels
addpath(genpath(path_Label));
%% user input file directory
path_Mask = uigetdir(pwd,'Select directory of 5D inferecned results'); % the directory where saves .nii 5D inferecned results
addpath(genpath(path_Mask));
SLabel = dir(fullfile(path_Label,'*.tif')); % pattern to match ground truth label filenames.
% unzip the files first
SMask = dir(fullfile(path_Mask,'*.nii')); % pattern to match output filenames.
if numel(SLabel) ~= numel(SMask)
    error('files numbers dont match.');
end
%% user input name of the saved file and directory 
prompt = {'Saved name', 'Output directory', 'confidence map: threshold', 'channel', 'dilation'};
dlgtitle = 'Input';
dims = [1 60];
definput = {'Evaluation_metrics.mat', 'seg_output/folder_name/','0.94', '2', '2' };
answer = inputdlg(prompt,dlgtitle,dims,definput);
eval_name = answer{1};
output_dir = answer{2};
thresh = str2double(answer{3});
channel = str2double(answer{4});
dilation = str2double(answer{5});
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
    figure_name = replace(SLabel(k).name, '_deconvolvedwholeexp_Label.tif', '');
    %% find the string name of SBR and density for plotting
    index = strfind(figure_name, '_');
    SBR_name(k,1) = str2double(figure_name(index(2)+1: index(3)-1));
    Density_name(k,1) =str2double( figure_name(index(4)+1: index(5)-1));
    figure_names(k,1) = str2double(figure_name(index(1)-1));
    %% start processing 
    GT_img = loadtiff(FLabel);
    FMask = fullfile(path_Mask,SMask(k).name);
    
    % third parameter determines the cutoff of the confidance map
    BW=process5Dimages_multiclass(FMask, thresh, channel);%0.94
    
    % second parameter determines the degree of img dilation, optional
    instance_seg_img = instance_segmentation_function(BW, dilation);
    %% Extract segmented cell information
    seg_cell_stats = regionprops3(instance_seg_img,'all');
    stats{k,1} = seg_cell_stats;
    %% generate a nii file with one integer indicating one instance fro testing purpose
    % comment out when not used
    try
        segName = strcat(output_dir,figure_name(1:end-9),'_seg.tif');
        write3Dtiff_V2(instance_seg_img, segName);
    catch
        mkdir(output_dir);
        segName = strcat(output_dir,figure_name(1:end-9),'_seg.tif');
        write3Dtiff_V2(instance_seg_img, segName);
    end
    %%
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
    %% Distance matching to calculate Distance of major axis, and relative angle of major axis.
    [distance_distribution,angle_distribution] = shortest_distance_matching_v2(GT_img, instance_seg_img);
    distance_distribution_eachframe{k,1} = distance_distribution;
    angle_distribution_eachframe{k,1} = angle_distribution;
    
    %%
    %compute the density for each GT image
    density(k,1) =  calculate_density(GT_img);
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
    density, FN_eachFrame, jaccard_counting_accuracy, FP_eachFrame,distance_distribution_eachframe, ...
    angle_distribution_eachframe);
% save the table as .mat in the matFileOutput folder
try
    full_eval_name = strcat('matFileOutput/',eval_name );
    save(full_eval_name, 'Evaluation_metrics');
catch 
    mkdir matFileOutput
    full_eval_name = strcat('matFileOutput/',eval_name );
    save(full_eval_name, 'Evaluation_metrics');
end

