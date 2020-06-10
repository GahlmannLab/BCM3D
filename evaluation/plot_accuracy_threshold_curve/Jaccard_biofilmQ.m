%% compute Jaccard for BiofilmQ results
% choose GT and 2T slices separately and compute accuracies based on
% different jaccrad compariosn threshold.
% choose Ground_truth Image
clear
[gt_files, gt_path] = uigetfile({'*.tif'},...
    'Open Images of ground truth',...
    'MultiSelect', 'on');

% choose Segment Image
[biofilmQresults_files, biofilmQresults_path] = uigetfile({'*.tif'},...
    'Open Images of Segment',...
    'MultiSelect', 'on');

jaccardRange = 0.1:0.1:0.9;
Jaccard_average = zeros(length(biofilmQresults_files), length(jaccardRange));
Dice_average = zeros(length(biofilmQresults_files), length(jaccardRange));
match_num = zeros(length(biofilmQresults_files), length(jaccardRange));
jaccard_countingAccuracy = zeros(length(biofilmQresults_files), length(jaccardRange));
dice_countingAccuracy = zeros(length(biofilmQresults_files), length(jaccardRange));
false_negative = zeros(length(biofilmQresults_files), length(jaccardRange));
false_positive = zeros(length(biofilmQresults_files), length(jaccardRange));

% h1 = figure;
% h2 = figure;

for j=1:length(biofilmQresults_files)
    j
    GT_img_path_temp=fullfile(gt_path, gt_files{j}); %get input_tiff with absolute path
    GT_img=loadtiff(GT_img_path_temp);
    
    biofilmQresult_path_temp=fullfile(biofilmQresults_path,biofilmQresults_files{j}); %get input_tiff with absolute path
    seg_img=loadtiff(biofilmQresult_path_temp);
    
%     size_data = size(GT_img);
%     for k=1:size_data(3)
%         figure(h1);
%         imagesc(GT_img(:,:,k));
%         figure(h2);
%         imagesc(seg_img(:,:,k));
%         pause
%     end

%     seg_img=process5Dimages(seg_img, 0.94); %for unlabelled data
%     seg_img = instance_segmentation_function(seg_img,2); %for unlabelled data
    
    % counting accuracy here is computed as jaccrard(TP/(TP+FP+FN))
    for i = 1:length(jaccardRange)
        [Jaccard_average(j, i), Dice_average(j, i), match_num(j, i), jaccard_countingAccuracy(j, i), ...
            dice_countingAccuracy(j, i), false_negative(j, i), false_positive(j, i)] = ...
            jaccard_Index_3D_Evaluation_v4_deleteMatched_JaccardAccuracy(GT_img, ...
            seg_img, jaccardRange(i));
    end
end

jaccard_countingAccuracy_avg = mean(jaccard_countingAccuracy);
dice_countingAccuracy_avg = mean(dice_countingAccuracy);

figure; hold on;
plot(jaccardRange, jaccard_countingAccuracy_avg);
xlabel('Jaccard index threshold');ylabel('Cell counting accuracy');
% xlabel('IoU threshold tau');ylabel('counting accuracy(Jaccard)');
% 
% plot(jaccardRange, dice_countingAccuracy_avg);
% legend('Jaccard Counting Accuracy','Dice Counting Accuracy');

% save('biofilmQ_500_1.mat')
% title('biofilmQ 500 1')