figure; hold on; box on;
all_centroid = sub_centroid;
all_contour = sub_contour;
for i = 1:length(all_centroid)
    color_ind = rand(1,3);
    color_ind = [0 0 1];

    plot3(all_centroid{i}(:,2),all_centroid{i}(:,1),all_centroid{i}(:,3),'.','color','r','MarkerSize',12);
    plot3(all_contour{i}(:,2),all_contour{i}(:,1),all_contour{i}(:,3),'-','color','b','LineWidth',1.5);
end
axis equal;
