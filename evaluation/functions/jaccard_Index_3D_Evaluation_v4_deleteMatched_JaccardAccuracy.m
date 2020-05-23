function [Jaccard_average, Dice_average, match_num, jaccard_countingAccuracy, dice_countingAccuracy, false_negative, false_positive] = jaccard_Index_3D_Evaluation_v4_deleteMatched_JaccardAccuracy(GT_img, seg_img, jaccardThresh)
    %%%segmentation evaluation used in instance segmentation
    %Jaccard index in 3D
    %parallel working increase speed
    %segmented objects have unique positive labels 
    %that are not necessarily propagated over time, 
    %background has zero label
    %Input: two 3D images; output: Jaccard_average(a number from 0 to 1)
    %2019-07-25
    %%
    %%set up parameters and maps
    l_gt = unique(GT_img);
    l_gt = l_gt(end);
    l_seg = unique(seg_img);
    num_seg = length(l_seg)-1;
    l_seg = l_seg(end);
    % set up parameter
    Jaccard_similarity = zeros(1,l_gt);
    Dice = zeros(1,l_gt);
    gt_index = cell(1,l_gt);
    for i=1:l_gt
    r1 = find(GT_img ==i);
    gt_index{i} = r1; 
    end
    seg_index = cell(1,l_seg);
    for j=1:l_seg
    r2 = find(seg_img == j);
    seg_index{j} = r2; 
    end
    %%
    %Match the segmented cell with GT
    % Old : try to match: criteria if they have 50% overlap thay are considered matched
    % Now: only match the ones with >0.5 jaccard, this prevents matching unsegmented
    % clusters. The unmatched will get a Jaccard of zero.
    for i = 1:l_gt
        for j = 1:l_seg
            inter = length(intersect(gt_index{1,i}, seg_index{1,j}));
            unionvar = length(union(gt_index{1,i}, seg_index{1,j}));
            if (inter/unionvar) > jaccardThresh
               %inter > (0.5 * length(gt_index{1,i}))
               Jaccard_similarity(1,i) = inter / unionvar;
               Dice(1,i) = 2 * Jaccard_similarity(1,i)/(Jaccard_similarity(1,i) + 1);
               gt_index{1,i} = NaN;
               seg_index{1,j} = NaN;
               break;
            end
        end
    end

    match_array = Jaccard_similarity(1,:);
    match_num = length(find(match_array>0));
    false_negative = length(gt_index) - match_num;
    %false_positive = length(seg_index)- match_num;
    false_positive = num_seg- match_num;
    % counting accuracy is jaccard not Dice anymore
    jaccard_countingAccuracy = (match_num) / (match_num + false_negative + false_positive);
    dice_countingAccuracy = (2*match_num) / (2*match_num + false_negative + false_positive);
    Jaccard_average = sum(match_array)/length(match_array);
    Dice_average = sum(Dice(1,:))/length(match_array);
end