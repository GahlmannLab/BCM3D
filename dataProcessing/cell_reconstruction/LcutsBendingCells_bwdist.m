% Mingxing Zhang, Gahlmann lab, Chemistry department, University of Virginia

function LcutsBendingCells_bwdist(file_path, file_name)
% load data
% V = loadtiff('normed_r1_Iter_12_150_2801_cropped__niftynet_out_segmented_0.95thresh_dilate2.tif');
% G = imread('normed_r1_Iter_12_150_2801_cropped__niftynet_out_segmented_0.95thresh_dilate2.tif', 1) ; % read in first image
% tiff_info = imfinfo('normed_r1_Iter_12_150_2801_cropped__niftynet_out_segmented_0.95thresh_dilate2.tif');
% %concatenate each successive tiff to tiff_stack
% for ii = 2 : size(tiff_info, 1)
%     temp_tiff = imread('normed_r1_Iter_12_150_2801_cropped__niftynet_out_segmented_0.95thresh_dilate2.tif', ii);
%     G = cat(3 , G, temp_tiff);
% end

% select under-segmented cells
% V_bad = V;
load (file_path);
beforeFit_V = V_noneedforcut;
beforeFit_clusterV = need_postV;
V_size = size(beforeFit_V);

% get large clusters' index
large_cluster_list = find(beforeFit_clusterV > 0);
[large_cluster_row, large_cluster_col, large_cluster_page] = ind2sub(size(beforeFit_clusterV), large_cluster_list);
large_cluster_list = [large_cluster_row, large_cluster_col, large_cluster_page];
num_good_cells = max(max(max(beforeFit_V)));

%% use the central axis obtained with LCuts to generate cells
post_post_segments = cell(length(post_segments), 1);
for k = 1:length(post_segments)
%     cell_axis_temp = post_segments{k};
%     bwdist_radius = post_radii{k};

    % fix bad axis (wrong order and duplicates)
    original_cell_axis_radius = [post_segments{k} post_radii{k}];
    processed_cell_axis_radius = fix_bad_cell_axis(original_cell_axis_radius);
    cell_axis_temp = double(processed_cell_axis_radius(:, 1:3));
    bwdist_radius = double(processed_cell_axis_radius(:, 4));
    
    % seems the radius are all smaller than the real data, so add 2 pixels
    bwdist_radius = bwdist_radius + 1.5;
    bwdist_index = post_radiidx{k};
    
    % figure; axis equal;
    % hold on;
    % for i = 1:length(cell_axis_temp)
    %     scatter3(cell_axis_temp(i, 1),cell_axis_temp(i, 2),cell_axis_temp(i, 3), 'r*');
    %     text(cell_axis_temp(i, 1),cell_axis_temp(i, 2),cell_axis_temp(i, 3), num2str(i));
    % end
    
    % Generate the normal vector for each circle
    cell_axis_normal = diff(cell_axis_temp);
    
    if length(cell_axis_normal) < 2
        points_cell_body = cell(2, 1);
    else
        points_cell_body = cell(length(cell_axis_normal) - 1, 1);
    end
    
    bin_size_circle = 0.01;
    theta_body = 0 : bin_size_circle : 2 * pi;
    
    % generate circles for cell body
    % radius = 4;
    % use the start position of each vector as the center
    [num_vectors, useless1] = size(cell_axis_normal);
    [num_points_onAxis, useless2] = size(cell_axis_temp);
    
    if num_vectors < 2
        for j = 1:num_points_onAxis
            center_body = cell_axis_temp(j, :);
            radius = bwdist_radius(j);
            normal_body = cell_axis_normal(1, :);
            v_body = null(normal_body);
            points_body = repmat(center_body', 1, size(theta_body,2)) + radius*(v_body(:,1)*cos(theta_body) + v_body(:,2)*sin(theta_body));
            points_cell_body{j} = points_body;
        end
    else
        for j = 2:num_vectors
            center_body = cell_axis_temp(j, :);
            radius = bwdist_radius(j);
            normal_body = cell_axis_normal(j - 1, :);
            v_body = null(normal_body);
            points_body = repmat(center_body', 1, size(theta_body,2)) + radius*(v_body(:,1)*cos(theta_body) + v_body(:,2)*sin(theta_body));
            points_cell_body{j - 1} = points_body;
            %     plot3(points_body(1,:), points_body(2,:), points_body(3,:),'r-');
        end
    end
    
    % % add bwdist result
    % data_size = size(post_V);
    % points_cell_body_bwdist = cell(length(points_cell_body) * 2, 1);
    % for i = 1:length(points_cell_body_bwdist)
    %     if rem(i, 2) == 0
    %         points_cell_body_bwdist{i} = points_cell_body{i/2};
    %     else
    %         bwdist_result = bwdist_index((i + 1) / 2);
    % %         bwdist_result_ind = zeros(length(bwdist_result), 3);
    %         [rows, col, pages] = ind2sub(data_size,bwdist_result);
    %         bwdist_result = [rows; col; pages];
    %         points_cell_body_bwdist{i} = bwdist_result;
    %     end
    % end
    
    
    % generate points for cell pole. The length of each cell pole is set as 4
    % pixels (the same as the radius of the cell). 4 even-distanced circles are drawn for the
    % cell pole. Their radius are 1/4, 1/2, 3/4 and 1 times the cell radius.
    num_circles_cell_pole = 5;
    normal_cell_pole = [0, 0, 1];
    cell_pole_points_temp = cell(num_circles_cell_pole, 1);
    cell_pole_1_points = cell(num_circles_cell_pole, 1);
    cell_pole_2_points = cell(num_circles_cell_pole, 1);
    
    % cell pole 1
    radius_pole1 = bwdist_radius(1);
    for i = 1:num_circles_cell_pole
        if i == 1
            angles = theta_body;
            center_temp = [0, 0, num_circles_cell_pole];
            radius_temp = 4;
            v_cell_pole = null(normal_cell_pole);
            points_temp = repmat(center_temp', 1, size(angles,2)) + radius_temp*(v_cell_pole(:,1)*cos(angles) + v_cell_pole(:,2)*sin(angles));
            %         points_temp = [0, 0, size_cell_pole];
        else
            bin_size_circle_pole = bin_size_circle * (i - 1);
            angles = 0:bin_size_circle_pole:2*pi;
            center_temp = [0, 0, ((num_circles_cell_pole - i) / (num_circles_cell_pole - 1)) * radius_pole1];
            radius_temp = ((i - 1) / (num_circles_cell_pole - 1)) * radius_pole1;
            v_cell_pole = null(normal_cell_pole);
            points_temp = repmat(center_temp', 1, size(angles,2)) + radius_temp*(v_cell_pole(:,1)*cos(angles) + v_cell_pole(:,2)*sin(angles));
        end
        cell_pole_1_points{i} = points_temp;
    end
    
    % cell pole 1
    radius_pole2 = bwdist_radius(end);
    for i = 1:num_circles_cell_pole
        if i == 1
            angles = theta_body;
            center_temp = [0, 0, num_circles_cell_pole];
            radius_temp = 4;
            v_cell_pole = null(normal_cell_pole);
            points_temp = repmat(center_temp', 1, size(angles,2)) + radius_temp*(v_cell_pole(:,1)*cos(angles) + v_cell_pole(:,2)*sin(angles));
            %         points_temp = [0, 0, size_cell_pole];
        else
            bin_size_circle_pole = bin_size_circle * (i - 1);
            angles = 0:bin_size_circle_pole:2*pi;
            center_temp = [0, 0, ((num_circles_cell_pole - i) / (num_circles_cell_pole - 1)) * radius_pole2];
            radius_temp = ((i - 1) / (num_circles_cell_pole - 1)) * radius_pole2;
            v_cell_pole = null(normal_cell_pole);
            points_temp = repmat(center_temp', 1, size(angles,2)) + radius_temp*(v_cell_pole(:,1)*cos(angles) + v_cell_pole(:,2)*sin(angles));
        end
        cell_pole_2_points{i} = points_temp;
    end
    
    % allocate space for the whole cell
    % whole_cell_points = cell(size_cell_pole * 2 + length(points_cell_body_bwdist), 1);
    whole_cell_points = cell(num_circles_cell_pole * 2 + length(points_cell_body), 1);
    
    % connect cell poles to cell body by rotation and moving
    original_vector_pole1 = normal_cell_pole;
    target_vector_pole1 = cell_axis_temp(1, :) - cell_axis_temp(2, :);
    original_vector_pole2 = normal_cell_pole;
    target_vector_pole2 = cell_axis_temp(end, :) - cell_axis_temp(end - 1, :);
    % rotate by vectors
    for i = 1 : length(cell_pole_1_points)
        % cell pole1
        cell_pole1_points_single = cell_pole_1_points{i};
        cell_pole1_points_single = rotCellWithVector(cell_pole1_points_single', original_vector_pole1, target_vector_pole1);
        cell_pole1_points_single = cell_pole1_points_single + cell_axis_temp(1, :);
        cell_pole_1_points{i} = cell_pole1_points_single';
        
        % cell pole2
        cell_pole2_points_single = cell_pole_2_points{i};
        cell_pole2_points_single = rotCellWithVector(cell_pole2_points_single', original_vector_pole2, target_vector_pole2);
        cell_pole2_points_single = cell_pole2_points_single + cell_axis_temp(end, :);
        cell_pole_2_points{i} = cell_pole2_points_single';
    end
    
    cell_pole_vertex = cell_pole_1_points{1};
    cell_pole_vertex = mean(cell_pole_vertex, 2);
    
    cell_pole_1_points{1} = cell_pole_vertex;
    % cell_pole_1_points{1} = [];
    
    cell_pole_vertex = cell_pole_2_points{1};
    cell_pole_vertex = mean(cell_pole_vertex, 2);
    
    cell_pole_2_points{1} =cell_pole_vertex;
    % cell_pole_2_points{1} = [];
    
    whole_cell_points(1:length(cell_pole_1_points)) = cell_pole_1_points(:);
    % whole_cell_points(length(cell_pole_1_points) + 1 : length(cell_pole_1_points) + length(points_cell_body_bwdist)) = points_cell_body_bwdist(:);
    whole_cell_points(length(cell_pole_1_points) + 1 : length(cell_pole_1_points) + length(points_cell_body)) = points_cell_body(:);
    cell_pole_2_points = flip(cell_pole_2_points);
    whole_cell_points(end - length(cell_pole_2_points) + 1 : end) = cell_pole_2_points;
    whole_cell_points_mat = [];
    
    % plot the whole cell
    % figure; axis equal; hold on;
    for j = 1:length(whole_cell_points)
        whole_cell_circles = whole_cell_points{j};
        %     if ~isempty(whole_cell_circles)
        %         plot3(whole_cell_circles(1,:), whole_cell_circles(2,:), whole_cell_circles(3,:),'b-');
        %     end
        whole_cell_points_mat = [whole_cell_points_mat; whole_cell_circles'];
    end
    
    % plot the cell axis
    % plot3(cell_axis_temp(:,1), cell_axis_temp(:,2), cell_axis_temp(:,3),'k-*');
    
    % apply these circles to extract cells from CNN confidence map
    whole_cell_points_mat = double(whole_cell_points_mat);
    in_cell_circles = inhull(large_cluster_list, whole_cell_points_mat);
    in_cell_circles_index = find(in_cell_circles);
    select_xyz = large_cluster_list(in_cell_circles_index,:);
    post_post_segments{k} = select_xyz;
    
    %     if length(select_xyz) < 100
    %         length(select_xyz)
    %         break;
    %     end
    % add post post processed cells into the good cell data
    for kk = 1:length(select_xyz)
        temp_index1 = select_xyz(kk, :);
        temp_index2 = sub2ind(V_size, temp_index1(1), temp_index1(2), temp_index1(3));
        beforeFit_V(temp_index2) = num_good_cells + k;
    end
    
    %     %     beforeFit_V(select_xyz) = num_good_cells + k;
    %     f1 = figure;
    %     hold on; axis equal;
    %     plot3(large_cluster_list(:,1), large_cluster_list(:,2), large_cluster_list(:,3),'r.');
    %     figure;hold on; axis equal;
    %     plot3(large_cluster_list(:,1), large_cluster_list(:,2), large_cluster_list(:,3),'r.');
    %     plot3(cell_axis_temp(:,1), cell_axis_temp(:,2), cell_axis_temp(:,3),'g-*');
    %     plot3(whole_cell_points_mat(:,1), whole_cell_points_mat(:,2), whole_cell_points_mat(:,3),'bo');
    %     pause;
    %     close(f1);
end

% relabel and filter small cells
num_cells = max(max(max(beforeFit_V)));
relabel_ind = 1;
for i = 1:num_cells
    if sum(sum(sum(beforeFit_V == i))) > 4 / 3 * pi * 64
        beforeFit_V(beforeFit_V == i) = relabel_ind;
        relabel_ind = relabel_ind + 1;        
    else
        beforeFit_V(beforeFit_V == i) = 0;
    end
end

post_post_seg_mat = beforeFit_V;
file_name_tif = strcat(file_name, '.tif');
save([file_name '.mat'], 'post_post_segments','post_post_seg_mat');
write3Dtiff_V2(uint16(beforeFit_V), file_name_tif);

% num_good_cells = max(max(max(beforeFit_V)));
% V_size = size(beforeFit_V);
% for kk = 1:length(post_post_segments)
%     select_xyz = post_post_segments{kk};
%     for i = 1:length(select_xyz)
%         temp_index1 = select_xyz(i, :);
%         temp_index2 = sub2ind(V_size, temp_index1(1), temp_index1(2), temp_index1(3));
%         beforeFit_V(temp_index2) = num_good_cells + kk;
%     end
% end

% bwdist;
end



