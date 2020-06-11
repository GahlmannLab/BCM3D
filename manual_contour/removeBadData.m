%% Remove bad contours and centroids
close all;
clear;

%% Load data and plot
[fileName,dirName] = uigetfile('*','Load the data');
dataFile = [dirName fileName];
cd (dirName);
data = load(dataFile);
folderName = fileName(1:end - 4);
mkdir(folderName);
cd (folderName);

cluster_centroid = data.group_1_centroid_cell;
cluster_contour = data.group_1_contour_cell;

f1 = figure; hold on;
for i = 1:length(cluster_centroid)
    centroid_1_temp = cluster_centroid{i};
    contour_1_temp = cluster_contour{i};
    
    color  = rand(1,3);
    
    centroid_1_temp_mat = centroid_1_temp(:,1:3);
    %     centroid_1_temp_mat = centroid_1_temp_mat(:,1:3);
    
    contour_1_temp_mat = contour_1_temp;
    
    % plot centroid and contours
    plot3(centroid_1_temp_mat(1),centroid_1_temp_mat(2),centroid_1_temp_mat(3),'o','color',color);
    text(centroid_1_temp_mat(1),centroid_1_temp_mat(2),centroid_1_temp_mat(3),num2str(i),'Color',color,'FontSize',14);
    
    plot3(contour_1_temp_mat(:,1),contour_1_temp_mat(:,2),contour_1_temp_mat(:,3),'s-','color',color);
    text(contour_1_temp_mat(1,1),contour_1_temp_mat(1,2),contour_1_temp_mat(1,3),num2str(i),'Color',color,'FontSize',14);
    
end

pause;

%% Determine if there are bad contours to be removed
answer1 = questdlg('Are there bad contours to remove?','Yes','No');
if strcmp(answer1,'Yes')
    
    answer2 = 'No';
    while strcmp(answer2,'No')
        cluster_centroid_temp = cluster_centroid;
        cluster_contour_temp = cluster_contour;
        % Select the centroid and contour to remove
        prompt = {'Select the index of the bad centroid'};
        title1 = 'Select bad centroids';
        dims = [2 50];
        %     definput = {'1'};
        in_para = inputdlg(prompt,title1,dims);
        index_bad_centroid = str2num(in_para{:});
        
        for j = 1:length(index_bad_centroid)
            cluster_centroid_temp{index_bad_centroid(j)} = [];
            cluster_contour_temp{index_bad_centroid(j)} = [];
        end
        
        f2 = figure; hold on;
        for i = 1:length(cluster_centroid_temp)
            if ~isempty(cluster_centroid_temp{i})
                
                centroid_1_temp2 = cluster_centroid_temp{i};
                contour_1_temp2 = cluster_contour_temp{i};
                
                color  = rand(1,3);
                
                centroid_1_temp_mat2 = centroid_1_temp2(:,1:3);
                %     centroid_1_temp_mat = centroid_1_temp_mat(:,1:3);
                
                contour_1_temp_mat2 = contour_1_temp2;
                
                % plot centroid and contours
                plot3(centroid_1_temp_mat2(1),centroid_1_temp_mat2(2),centroid_1_temp_mat2(3),'o','color',color);
                text(centroid_1_temp_mat2(1),centroid_1_temp_mat2(2),centroid_1_temp_mat2(3),num2str(i),'Color',color,'FontSize',14);
                
                plot3(contour_1_temp_mat2(:,1),contour_1_temp_mat2(:,2),contour_1_temp_mat2(:,3),'s-','color',color);
                text(contour_1_temp_mat2(1,1),contour_1_temp_mat2(1,2),contour_1_temp_mat2(1,3),num2str(i),'Color',color,'FontSize',14);
            else
%                 close (f2);
                continue;
            end
        end
        
        pause;
        answer2 = questdlg('Do you satisfy with this good cluster?','Yes','No');
        
        if strcmp(answer2,'Yes')
            cluster_centroid_temp = cluster_centroid_temp(~cellfun('isempty',cluster_centroid_temp));
            group_1_centroid_cell = cluster_centroid_temp;
            
            cluster_centroid_temp_mat = cell2mat(cluster_centroid_temp);
            cluster_centroid_temp_mat = cluster_centroid_temp_mat(:,1:3);
            good_cluster_temp = cluster_centroid_temp_mat;
            
            cluster_contour_temp = cluster_contour_temp(~cellfun('isempty',cluster_contour_temp));
            group_1_contour_cell = cluster_contour_temp;
            
            good_contour_temp_mat = cell2mat(cluster_contour_temp);
            good_contour_temp = good_contour_temp_mat;
            
            save([folderName '_manual_removed' '.mat'], 'good_cluster_temp','good_contour_temp','group_1_contour_cell','group_1_centroid_cell');
            
        else
            continue;
        end
    end
    close all;
    clear;
else
    disp('Good data, no bad contours.');
    close all;
    clear;
end
