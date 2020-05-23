% This file includes two parts:
% Part 1: visualize the data in matlab
% Part 2: compare the result with ground truth use shortest distance matching.
% Input: ground truth centroid: myx_all_centroid
%        LCUTS result: final_Comp

% Last Edit: Jie Wang, 2019-05-30
%%
clc;clear all;close all;


%% load data
load fitted_cell_points.mat;


load synthetic_rotation_shewan.mat;
% myx_all_centroid
% swap axis if needed
final_Comp = new_cell.';

%% ------------------Part 1 data visualization----------------------------
figure;
hold on;axis equal;
cmap = jet(size(final_Comp,2));
for q = 1:size(final_Comp,2)
    colorm = rand(1,3);
    y = final_Comp{1,q}(:,2);
    x = final_Comp{1,q}(:,1);
    z = final_Comp{1,q}(:,3);
    
    %plot3(x,y,z,'.');
    plot3(x,y,z,'.','Color',colorm,'MarkerSize',10)
    axis([0 300 0 300 ]);
    %plot(x,y,'.');
    %pause;
    hold on;
end

%% reshape the groundtruth
stat_gt = zeros(size(fitted_model,1),6);
countnum = 1;
figure;
axis([0 300 0 300 ])
for i = 1: size(fitted_model,1)
   currentSeg = fitted_model{i,1}.';
   currentSeg(currentSeg(:,1)>300 | currentSeg(:,2)>300 | currentSeg(:,3)<20 | currentSeg(:,3)>91,:) = [];
   currentSeg(:,3) = currentSeg(:,3)-20;
   if ~isempty(currentSeg)
       hold on;plot3(currentSeg(:,1),currentSeg(:,2),currentSeg(:,3),'.');
       maxrange = max(currentSeg,[],1);
       minrange = min(currentSeg,[],1);
       center = mean([maxrange;minrange],1);
       [orientation] = GetCompOrientation(currentSeg);
       stat_gt(i,1:3) = center;
       stat_gt(i,4:6) = orientation;
       gt_fixedpoints{countnum} =  currentSeg;
       countnum = countnum+1;
   else
       stat_gt(i,1:6) = -20;
   end
end
stat_gt(stat_gt==-20) = [];
stat_gt = reshape(stat_gt,[352,6]);



%% ------------------Part 2 comparison ------------------------------------
%% statistics from ground truth
%stat_gt = zeros(size(manual_centroid,2),6);
%for i = 1: size(manual_centroid,1)
    %currentSeg = manual_centroid{i,1};
    %center = mean(currentSeg,1);
    %orientation = GetCompOrientation(currentSeg);
    %stat_gt(i,1:3) = center;
    %stat_gt(i,4:6) = orientation;
%end

%% statistics from experiment data
stat_e = zeros(size(final_Comp,2),6);
for j = 1: size(final_Comp,2)
    currentSeg = final_Comp{1,j};
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
%% Redo former two steps if needed.


%% show: clean-up
correspond(unmatched_inG)=0;
correspondence = zeros(length(final_Comp),2);
for num = 1:length(final_Comp)
     if correspond(num) == 0
        correspondence(num,1) = nan;
        correspondence(num,2) = nan;
    else
        correspondence(num,1) = dists(num,correspond(num));
        correspondence(num,2) = real(theta(num,correspond(num)));
    end
end
%%

% figure;histogram(correspondence(:,1),20);title('distance cost by best distance matching ');
% figure;histogram(correspondence(:,2),20,'BinLimits',[0,90]);title('relative angle by best distance matching');

%% show with markers
h1 = histogram(correspondence(:,1),20,'BinLimits',[0,20]);
counts1 = h1.Values;
E1 = h1.BinEdges;
xloc1 = E1(1:end-1)+diff(E1)/2;

figure;histogram(correspondence(:,1),20,'BinLimits',[0,20]);title('Distance cost by best distance matching');ylim([0 80]);
text(xloc1, counts1, num2str(counts1'),'HorizontalAlignment','center', 'VerticalAlignment','bottom');


figure;
h2 = histogram(correspondence(:,2),20,'BinLimits',[0,90]);
counts2 = h2.Values;
E2 = h2.BinEdges;
xloc2 = E2(1:end-1)+diff(E2)/2;

figure;histogram(correspondence(:,2),20,'BinLimits',[0,90]);title('Relative angle by best distance matching');ylim([0 130]);
text(xloc2, counts2, num2str(counts2'),'HorizontalAlignment','center', 'VerticalAlignment','bottom');


%% tools
function direction = GetCompOrientation(currentSeg)
% This function is to find the direction of current component
Sub = ones(size(currentSeg));
Sub = mean(currentSeg,1).*Sub;
ShiftSeg = currentSeg - Sub;

CovM = cov(ShiftSeg);
[V,D] = eig(CovM);
x_v1= V(1,3);
y_v1= V(2,3);
z_v1= V(3,3);
direction = [x_v1,y_v1,z_v1];
end

function y = inner(a,b)
% This is a MatLab function to compute the inner product of
% two vectors a and b.
% By Ralph Howard on 1/12/98
% Call syntax: y = inner(a,b) or inner(a,b)
% Input: The two vectors a and b
% Output: The value of the inner product of a and b.
c=0;            % intialize the variable c
n= length(a);        % get the lenght of the vector a
for k=1:n        % start the loop
    c=c+a(k)*b(k);    % update c by the k-th product in inner product
end
y = c;
end