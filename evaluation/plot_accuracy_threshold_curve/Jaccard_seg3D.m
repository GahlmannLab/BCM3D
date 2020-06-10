%% compute Jaccard for BiofilmQ results
% choose GT and 2T slices separately and compute accuracies based on
% different jaccrad compariosn threshold.
% choose Ground_truth Image
clear
[gt_files, gt_path] = uigetfile({'*.tif'},...
    'Open Images of ground truth',...
    'MultiSelect', 'on');

% choose Segment Image
[biofilmQresults_files, biofilmQresults_path] = uigetfile({'*.mat'},...
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
    
    %     GT_img = open_3D_tiff(GT_img_path_temp);
    %     GT_img=loadtiff(GT_img_path_temp);
    
    %     GT_img=load(GT_img_path_temp);
    
    GT_img = imread(GT_img_path_temp, 1) ; % read in first image
    tiff_info = imfinfo(GT_img_path_temp);
    %concatenate each successive tiff to tiff_stack
    for ii = 2 : size(tiff_info, 1)
        temp_tiff = imread(GT_img_path_temp, ii);
        GT_img = cat(3 , GT_img, temp_tiff);
    end
    
    biofilmQresult_path_temp=fullfile(biofilmQresults_path,biofilmQresults_files{j}); %get input_tiff with absolute path
    load(biofilmQresult_path_temp);
    seg_img = total_cells;
    
    %     seg_img = niftiread(biofilmQresult_path_temp);
    %     seg_img = process5Dimages(biofilmQresult_path_temp, 0.94);
    %     seg_img = instance_segmentation_function(seg_img, 2); %for unlabelled data
    %     seg_img = imrotate(uint16(seg_img), 90);
    %     seg_img = flip(seg_img, 1);
    
%     size_data = size(GT_img);
%     for k=1:size_data(3)
%         figure(h1);
%         imagesc(GT_img(:,:,k));
%         figure(h2);
%         imagesc(seg_img(:,:,k));
%         pause
%     end
%     
%     size_data = size(GT_img);
%     for k=1:size_data(1)
%         figure(h1);
%         imagesc(squeeze(GT_img(k,:,:)));
%         figure(h2);
%         imagesc(squeeze(seg_img(k,:,:)));
%         pause
%     end
    
    % counting accuracy here is computed as jaccrard(TP/(TP+FP+FN))
    for i = 1:length(jaccardRange)
        [Jaccard_average(j, i), Dice_average(j, i), match_num(j, i), jaccard_countingAccuracy(j, i), ...
            dice_countingAccuracy(j, i), false_negative(j, i), false_positive(j, i)] = ...
            jaccard_Index_3D_Evaluation_v4_deleteMatched_JaccardAccuracy(GT_img, ...
            seg_img, jaccardRange(i));
    end
    
    
    %     for i = 1:length(jaccardRange)
    %         [match_num(j, i), false_positive(j, i), false_negative(j, i), Jaccard_average(j, i), Dice_average(j, i)] = countingAcc_test(GT_img, seg_img, jaccardRange(i));
    %Dice(i) =  2*TP(i)/(2*TP(i)+FP(i)+FN(i));
    %     end
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

% save('seg3D_2000_1.mat')
% title('seg3D 2000 1')