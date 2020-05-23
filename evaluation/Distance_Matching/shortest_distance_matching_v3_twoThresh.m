function [jaccard_countingAccuracy] = shortest_distance_matching_v3_twoThresh(GT_img, seg_img, thetaThresh, distThresh)

    %% ------------------Part 1 setting and change matrix to cell array----------------------------
    %%set up parameters and maps
    l_gt = unique(GT_img);
    l_gt = l_gt(end);
    l_seg = unique(seg_img);
    l_seg = l_seg(end);

    % set up parameter
    gt_index = cell(l_gt,1);
    for i=1:l_gt
    [x1,y1,z1] = ind2sub(size(GT_img),find(GT_img ==i));
    gt_index{i} = [x1,y1,z1]; 
    end
    seg_index = cell(l_seg,1);
    for j=1:l_seg
    [x2,y2,z2] = ind2sub(size(seg_img),find(seg_img ==j));
    seg_index{j} = [x2,y2,z2]; 
    end
    
    gt_index = gt_index(~cellfun(@isempty, gt_index(:,1)), :);
    seg_index = seg_index(~cellfun(@isempty, seg_index(:,1)), :);
 
    %% ------------------Part 2 comparison ------------------------------------
    %% statistics from ground truth
    stat_gt = zeros(size(gt_index,2),6);
    for i = 1: size(gt_index,1)
        currentSeg = gt_index{i,1};
        center = mean(currentSeg,1);
        orientation = GetCompOrientation(currentSeg);
        stat_gt(i,1:3) = center;
        stat_gt(i,4:6) = orientation;
    end

    %% statistics from experiment data
    stat_e = zeros(size(seg_index,2),6);
    for j = 1: size(seg_index,1)
        currentSeg = seg_index{j,1};
        center = mean(currentSeg,1);
        orientation = GetCompOrientation(currentSeg);
        stat_e(j,1:3) = center;
        stat_e(j,4:6) = orientation;
    end

    %% calculate relative angle
    dists = zeros(length(stat_gt(:,1)),length(stat_e(:,1)));
    theta = zeros(size(dists,1),size(dists,2)); 
    for o = 1:size(dists,1)
        for oo = 1: size(dists,2)
            a = stat_gt(o,4:6);
            b = stat_e(oo,4:6);
            innervalue = inner(a,b);
            costheta(o,oo) = abs(innervalue)/(norm(a,2)*norm(b,2));
            theta(o,oo) = 180*acos(costheta(o,oo))/pi;
        end
    end
    %% order centroids
    dists = zeros(length(stat_gt(:,1)),length(stat_e(:,1))); %according to the size of the component
    sdists = zeros(length(stat_gt(:,1)),length(stat_e(:,1)));
    map = zeros(length(stat_gt(:,1)),length(stat_e(:,1)));
    for num = 1:size(stat_gt,1)
        currentDist = stat_e(:,1:3) - stat_gt(num,1:3);
        edist = vecnorm(currentDist,2,2);
        %current_idx = find(edist == min(edist));
        [sedist,s_ind] = sort(edist);
        dists(num,:) = edist.';
        sdists(num,:) = sedist.';
        map(num,:) = s_ind.';   
    end
    
    %% threshold
    matchedCells = zeros(1,length(gt_index));
     for i = 1:size(dists,1)
        for j = 1: size(dists,2)
            if theta(i,j) < thetaThresh && sdists(i,j)< distThresh
               matchedCells(1,i) = 1;
               theta(i,j) = NaN;
               sdists(i,j) = NaN;
            end
        end
     end
     
    match_num = length(find(matchedCells>0));
    false_negative = size(dists,1) - match_num;
    %false_positive = length(seg_index)- match_num;
    false_positive = size(dists,2) - match_num;
    % counting accuracy is jaccard not Dice anymore
    jaccard_countingAccuracy = (match_num) / (match_num + false_negative + false_positive);
end