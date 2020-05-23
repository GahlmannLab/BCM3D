function Clustering_1st()
clear all
close all

[fileName,dirName] = uigetfile('*','Load the centoid data');
dataFile = [dirName fileName];
data = load(dataFile);
folderName = fileName(1:end - 4);
mkdir(folderName);
cd (folderName);

cluster_2_centroid = data.good_cluster_temp;
cluster_2_contour = data.good_contour_temp;
% x = [];
% y = [];
% z = [];
x = zeros(length(cluster_2_centroid),1);
y = zeros(length(cluster_2_centroid),1);
z = zeros(length(cluster_2_centroid),1);

% for i = 1:length(cluster_2_centroid)
%     x = [x;cluster_2_centroid{i}(1)];
%     y = [y;cluster_2_centroid{i}(2)];
%     z = [z;cluster_2_centroid{i}(3)];
% end

for i = 1:length(cluster_2_centroid)
    x(i) = cluster_2_centroid{i}(1);
    y(i) = cluster_2_centroid{i}(2);
    z(i) = cluster_2_centroid{i}(3);
end
X = [x,y,z];

% plot the original data, all centroids
scrsz = get(0,'ScreenSize');
f1 = figure; hold on;
set(f1,'Position',[scrsz(3)/2-300 scrsz(4)/2 550 450]);
scatter3(x,y,z,'r.');
axis square;

f2 = figure; hold on;
set(f2,'Position',[scrsz(3)/2+300 scrsz(4)/2 550 450]);
for m = 1:length(cluster_2_contour)
    % color_ind = rand(1,3);
    all_contours_f2 = plot3(cluster_2_contour{m}(:,1),cluster_2_contour{m}(:,2),cluster_2_contour{m}(:,3),'g.-');
    all_centroids_f2 = plot3(cluster_2_centroid{m}(1),cluster_2_centroid{m}(2),cluster_2_centroid{m}(3),'g.-');
end
centroid_clusters = {};
contour_clusters = {};
k = 1;
epsilon = 4;
MinPts = 6;
answer2 = 'Yes';
% determine which cluster to be checked in detail
while strcmp(answer2,'Yes')
    [IDX, isnoise, answer1, X_new,cluster_2_contour_new,cluster_2_centroid_new,good_cluster_temp,good_contour_temp] = DBSCAN_fine_sel(X,epsilon,MinPts,f1,f2,cluster_2_contour,cluster_2_centroid,k);
    save(['newX_' num2str(k) '.mat'], 'X_new');
    centroid_clusters{k} = good_cluster_temp;
    contour_clusters{k} = good_contour_temp;
    answer2 = questdlg('Do you want to continue clustering centroids?','Yes','No');
    X = X_new;
    cluster_2_contour = cluster_2_contour_new;
    cluster_2_centroid = cluster_2_centroid_new;
    k = k + 1;
end
save('semiautocluster.mat','centroid_clusters','contour_clusters');
end
