
% Part 1: compare the result with ground truth using shortest distance matching.
% Input: ground truth centroid: seg_img, 3D matrix
%        Ground truth result: GT_img, 3d matrix
% predictionmask: instance_seg_img
% modified based on Jie's code 2019-05-30 on 2019-07-30 Yibo
% 
%%
function [distance_distribution,angle_distribution] = shortest_distance_matching_v2(GT_img, seg_img)

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

    %%
    correspond = map(:,1);

    %% correct correspondence
    check_dists = sdists(:,1);
    %correspond(check_dists>20)=0;
    sum(check_dists>20)
    %% Find out duplicated matching
    clear unsure_idx;
    count = 1;
    needmatch_inE = [];
    for c = 1:size(stat_e,1)
        idx = find(correspond == c);
        if (size(idx,1)~=1) && (~isempty(idx))
            unsure_idx{count} = [idx;c];
            count = count + 1;
        elseif isempty(idx)
            needmatch_inE = [needmatch_inE;c];
        end
    end

    unmatched_inG = [];
    for co = 1: count -1
        currentidx = unsure_idx{1,co}(1:end-1);
        comparetheta = theta(currentidx,unsure_idx{1,co}(end));
        comparedist = dists(currentidx,unsure_idx{1,co}(end));
        comparevalue =  comparetheta.*comparedist;
        good_idx = find(comparevalue == min(comparevalue(:)));
        unmatched_idx = currentidx(currentidx~=currentidx(good_idx));
        unmatched_inG = [unmatched_inG;unmatched_idx];
    end

    %% try to match unmatched
    try
        unmatched_dists = dists(unmatched_inG,needmatch_inE); 
        % sort
        [s_unmatched_dists,sortE_map] = sort(unmatched_dists,2); % the order of E
        [shortest_matches,shortest_idx] = sort(s_unmatched_dists(:,1)); % the order of G

        smallnum = min([size(unmatched_inG,1),size(needmatch_inE,1)]);
        index_tomatch = sortE_map(shortest_idx(1:smallnum),1);
        % for g = 1: smallnum%size(unmatched_inG)
        %     correspond(unmatched_inG(index_tomatch(g))) = needmatch_inE(sortE_map(index_tomatch(g),1));
        % end
        for g = 1: smallnum%size(unmatched_inG)
            correspond(unmatched_inG(index_tomatch(g))) = needmatch_inE(sortE_map(index_tomatch(g),1));
        end
    catch 
        warning('all cells are matched');
    end
    %% Redo former two steps if needed.


    %% show: clean-up
    correspond(unmatched_inG)=0;
    correspondence = zeros(length(seg_index),2);  
    if length(seg_index) > length(gt_index)
        smaller_num = length(gt_index);
    else
        smaller_num = length(seg_index);
    end
    
    for num = 1:smaller_num
         if correspond(num) == 0
            correspondence(num,1) = nan;
            correspondence(num,2) = nan;
        else
            correspondence(num,1) = dists(num,correspond(num));
            correspondence(num,2) = real(theta(num,correspond(num)));
        end
    end
    

    %% Plot and show with markers

%     h1 = histogram(correspondence(:,1),20,'BinLimits',[0,20]);
%     counts1 = h1.Values;
%     E1 = h1.BinEdges;
%     xloc1 = E1(1:end-1)+diff(E1)/2;
     distance_distribution = correspondence(:,1);
%     figure;
%     histogram(correspondence(:,1),20,'BinLimits',[0,20]);
%     title('Best distance matching: Distance distribution compariosn');
%     ylim([0 80]);
%     text(xloc1, counts1, num2str(counts1'),'HorizontalAlignment','center', 'VerticalAlignment','bottom');
% 
     angle_distribution = correspondence(:,2);
% 
%     h2 = histogram(correspondence(:,2),20,'BinLimits',[0,90]);
%     counts2 = h2.Values;
%     E2 = h2.BinEdges;
%     xloc2 = E2(1:end-1)+diff(E2)/2;
% 
%     figure;
%     histogram(correspondence(:,2),20,'BinLimits',[0,90]);
%     title('Best distance matching: Relative angle distribution');ylim([0 130]);
%     text(xloc2, counts2, num2str(counts2'),'HorizontalAlignment','center', 'VerticalAlignment','bottom');
end

