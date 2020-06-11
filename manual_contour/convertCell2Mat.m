% Convert the cell format data to matrix

[fileName,dirName] = uigetfile('*','Load the centoid data');
dataFile = [dirName fileName];
data = load(dataFile);
folderName = fileName(1:end - 4);
mkdir(folderName);
cd (folderName);

all_cells_centroid = data.cells_centroid_good_nonEmpty;
all_cells_contour = data.cells_contours_good_nonEmpty;

% get the total number of elements in the centroid and contour data
length_all_centroid = 0;
length_all_contours = 0;

for j = 1:length(all_cells_centroid)
    
    length_all_centroid = length_all_centroid + length(all_cells_centroid{j});
    all_cells_contour_temp = all_cells_contour{j};
    length_all_contours_temp = 0;
    for k = 1:length(all_cells_contour_temp)
        length_all_contours_temp = length_all_contours_temp + length(all_cells_contour_temp{k});
    end
    length_all_contours = length_all_contours + length_all_contours_temp;
end

% allocate memory for all centroid and contours
all_centroid_mat = cell(0);
all_contours_mat = cell(0);

for i = 1:length(all_cells_centroid)
    
    cell_centroid_cel_temp = all_cells_centroid{i};
    cell_centroid_mat_temp = cell2mat(cell_centroid_cel_temp);
    cell_centroid_mat_temp = cell_centroid_mat_temp(:,1:3);
    for ii = 1:length(cell_centroid_cel_temp)
        temp_centroid = cell_centroid_mat_temp(ii,:);
        all_centroid_mat = [all_centroid_mat; temp_centroid];
    end
%     all_centroid_mat = [all_centroid_mat; cell_centroid_mat_temp];
    
    cell_contour_cel_temp = all_cells_contour{i};
%     cell_contour_mat_temp = cell2mat(all_cells_contour{i});
    
    for jj = 1:length(cell_contour_cel_temp)
        temp_contour = cell_contour_cel_temp{jj};
        all_contours_mat = [all_contours_mat; temp_contour];
    end    
%     all_contours_mat = [all_contours_mat; cell_contour_mat_temp]; 
    sub_centroid = cell_centroid_cel_temp;
    sub_contour = cell_contour_cel_temp;
    save([folderName ' subCell ' num2str(i) '.mat'], 'sub_centroid','sub_contour');
    
end
