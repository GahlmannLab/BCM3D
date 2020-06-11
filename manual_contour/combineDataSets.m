% Combine data sets
[fileName,dirName] = uigetfile('*','Load the centoid data','MultiSelect','on');

folderName = fileName{1};
folderName = folderName(1:end - 4);
folderName = [folderName '_combined'];
mkdir(folderName);
cd (folderName);
good_cluster_temp_temp = [];
good_contour_temp_temp = [];
group_1_centroid_cell_temp = cell(0);
group_1_contour_cell_temp = cell(0);

myx_all_centroid = cell(length(fileName),1);
myx_all_contour = cell(length(fileName),1);

for i = 1:length(fileName)
    dataFile = [dirName fileName{i}];
    data = load(dataFile);
    
    myx_all_centroid{i} = data.good_cluster_temp;
    myx_all_contour{i} = data.good_contour_temp;
    
    good_cluster_temp_temp = [good_cluster_temp_temp; data.good_cluster_temp];
    good_contour_temp_temp = [good_contour_temp_temp; data.good_contour_temp];
    
    group_1_centroid_cell_temp = [group_1_centroid_cell_temp; data.group_1_centroid_cell];
    group_1_contour_cell_temp = [group_1_contour_cell_temp; data.group_1_contour_cell];
    
end

good_cluster_temp = good_cluster_temp_temp;
good_contour_temp = good_contour_temp_temp;
group_1_centroid_cell = group_1_centroid_cell_temp;
group_1_contour_cell = group_1_contour_cell_temp;

save([folderName '.mat'],'good_cluster_temp','good_contour_temp','group_1_centroid_cell','group_1_contour_cell','myx_all_centroid','myx_all_contour');

