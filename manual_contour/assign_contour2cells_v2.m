% assign contour to cells, the idea is to calculate the controid of each
% contour in each direction, link these controid to generate three axis,
% find the cross point of these three axis, use this point and three axis
% to generate a ellipsoid, estimate if a contour belongs to this ellipsoid
figure;hold on;title('Original Contours'); 
x_coord = 1:1:length(x_frames);
contour_data = cell(length(x_frames) + length(y_frames) + length(z_frames),1);
x_data = cell(length(x_frames),1);
y_data = cell(length(y_frames),1);
z_data = cell(length(z_frames),1);
for i = 1:length(x_frames)
    d = 1; % means frames in x direction
    for j = 1:length(x_frames(i).cellContour)
        if ~isempty(x_frames(i).cellContour{j}) && length(x_frames(i).cellContour{j}) > 4
            y = x_frames(i).cellContour{j}(:,2); % the y coordinates in each frame is the y in the original data
            z = x_frames(i).cellContour{j}(:,1); % the x coordinates in each frame is the z in the original data
            x = x_coord(i)*ones(length(y),1);
            polyin = polyshape({z},{y});
            [zc,yc] = centroid(polyin);
            xc = x_coord(i);
            x_data{i}{j,1} = [x y z]; % save coordinates of each contour
            x_data{i}{j,2} = [xc yc zc d i j]; % save coordinates of the corresponding centroid (xc yc zc),
                                           % the direction (d = 1-x, 2-y, 3-z), frame number (i), contour number (j)            
%             plot(polyin)
            plot3(x,y,z,'r-*')
            plot3(xc,yc,zc,'ro')
%             clear y z;
        else
        end
    end
end

% y frames, the xz plane of the original data set, need to assign y values for each
% frame 
y_coord = 1:1:length(y_frames);
for i = 1:length(y_frames)
    d = 2; % means frames in y direction
    for j = 1:length(y_frames(i).cellContour)
        if ~isempty(y_frames(i).cellContour{j}) && length(y_frames(i).cellContour{j}) > 4
            z = y_frames(i).cellContour{j}(:,1); % the x coordinates in each frame is the z in the original data
            x = y_frames(i).cellContour{j}(:,2); % the y coordinates in each frame is the x in the original data
            y = y_coord(i)*ones(length(x),1);
            polyin = polyshape({z},{x});
            [zc,xc] = centroid(polyin);
            yc = y_coord(i);
            y_data{i}{j,1} = [x y z]; % save coordinates of each contour
            y_data{i}{j,2} = [xc yc zc d i j]; % save coordinates of the corresponding centroid (xc yc zc),
                                           % the direction (d = 1-x, 2-y, 3-z), frame number (i), contour number (j)  
%             plot(polyin)
            plot3(x,y,z,'b-*')
            plot3(xc,yc,zc,'bo')

%             clear y z;
        else           
        end
    end
end

% z frames, the xy plane of the original data set, need to assign z values for each
% frame 
z_coord = 1:1:length(z_frames);
for i = 1:length(z_frames)
    d = 3; % means frames in z direction
    for j = 1:length(z_frames(i).cellContour)
        if ~isempty(z_frames(i).cellContour{j}) && length(z_frames(i).cellContour{j}) > 4
            x = z_frames(i).cellContour{j}(:,1); % the x coordinates in each frame is the x in the original data
            y = z_frames(i).cellContour{j}(:,2); % the y coordinates in each frame is the y in the original data
            z = z_coord(i)*ones(length(x),1);
            polyin = polyshape({x},{y});
            [xc,yc] = centroid(polyin);
            zc = z_coord(i);
            z_data{i}{j,1} = [x y z]; % save coordinates of each contour
            z_data{i}{j,2} = [xc yc zc d i j]; % save coordinates of the corresponding centroid (xc yc zc),
                                           % the direction (d = 1-x, 2-y, 3-z), frame number (i), contour number (j)  
%             plot(polyin)
            plot3(x,y,z,'g-*')
            plot3(xc,yc,zc,'go')
%             clear y z;
        else
        end
    end
end
contour_data = [x_data; y_data; z_data];
% remove empty element
contour_data = contour_data(~cellfun('isempty',contour_data));
contour_data = cat(1,contour_data{:});
contour_data = contour_data(~cellfun('isempty',contour_data(:,1)),:); % remove empty element again
contourOnly_belong = cell(length(contour_data));
cell_inf = cell(length(contour_data),length(contour_data),2);
contour_matrix = zeros(length(contour_data));

for i = 1:length(contour_data)
    contour_temp_0 = contour_data{i,1}; % the centroid of the first contour
    centroid_temp_0 = contour_data{i,2}(1:3); % the centroid of the first contour
    frame_direction_0 = contour_data{i,2}(4);
    
    contour_belong = cell(length(contour_data),2);
    for j = 1:length(contour_data)
        switch frame_direction_0
            case 1 % contour_temp_0 is in the x direction
                if contour_data{j,2}(4) == 1 % project along x direction to yz plane
%                     in = inpolygon(centroid_temp_0(2),centroid_temp_0(3),contour_data{j,1}(:,2),contour_data{j,1}(:,3));
%                     % check back if the centroid of the contour found belongs to
%                     % the contour for the original centroid
%                     if in == 1
%                         centroid_temp_1 = contour_data{j,2}(1:3); % centroid of the found contour
%                         % project centroid_temp_1 along x direction to yz plane
%                         % where the contour_temp_0 stays.
%                         in2 = inpolygon(centroid_temp_1(2),centroid_temp_1(3),contour_temp_0(:,2),contour_temp_0(:,3));
%                     else
                        in2 = 0;
%                     end
                elseif contour_data{j,2}(4) == 2 % project along y direction to xz plane
                    in = inpolygon(centroid_temp_0(1),centroid_temp_0(3),contour_data{j,1}(:,1),contour_data{j,1}(:,3));
                    % check back if the centroid of the contour found belongs to
                    % the contour for the original centroid
                    if in == 1
                        centroid_temp_1 = contour_data{j,2}(1:3); % centroid of the found contour
                        % project centroid_temp_1 along x direction to yz plane
                        % where the contour_temp_0 stays.
                        in2 = inpolygon(centroid_temp_1(2),centroid_temp_1(3),contour_temp_0(:,2),contour_temp_0(:,3));
                    else
                        in2 = 0;
                    end
                elseif contour_data{j,2}(4) == 3 % project along z direction to xy plane
                    in = inpolygon(centroid_temp_0(1),centroid_temp_0(2),contour_data{j,1}(:,1),contour_data{j,1}(:,2));
                    % check back if the centroid of the contour found belongs to
                    % the contour for the original centroid
                    if in == 1
                        centroid_temp_1 = contour_data{j,2}(1:3); % centroid of the found contour
                        % project centroid_temp_1 along x direction to yz plane
                        % where the contour_temp_0 stays.
                        in2 = inpolygon(centroid_temp_1(2),centroid_temp_1(3),contour_temp_0(:,2),contour_temp_0(:,3));
                    else
                        in2 = 0;
                    end
                end
            case 2 % contour_temp_0 is in the y direction
                if contour_data{j,2}(4) == 1 % project along x direction to yz plane
                    in = inpolygon(centroid_temp_0(2),centroid_temp_0(3),contour_data{j,1}(:,2),contour_data{j,1}(:,3));
                    % check back if the centroid of the contour found belongs to
                    % the contour for the original centroid
                    if in == 1
                        centroid_temp_1 = contour_data{j,2}(1:3); % centroid of the found contour
                        % project centroid_temp_1 along y direction to xz plane
                        % where the contour_temp_0 stays.
                        in2 = inpolygon(centroid_temp_1(1),centroid_temp_1(3),contour_temp_0(:,1),contour_temp_0(:,3));
                    else
                        in2 = 0;
                    end
                elseif contour_data{j,2}(4) == 2 % project along y direction to xz plane
%                     in = inpolygon(centroid_temp_0(1),centroid_temp_0(3),contour_data{j,1}(:,1),contour_data{j,1}(:,3));
%                     % check back if the centroid of the contour found belongs to
%                     % the contour for the original centroid
%                     if in == 1
%                         centroid_temp_1 = contour_data{j,2}(1:3); % centroid of the found contour
%                         % project centroid_temp_1 along y direction to xz plane
%                         % where the contour_temp_0 stays.
%                         in2 = inpolygon(centroid_temp_1(1),centroid_temp_1(3),contour_temp_0(:,1),contour_temp_0(:,3));
%                     else
                        in2 = 0;
%                     end
                elseif contour_data{j,2}(4) == 3 % project along z direction to xy plane
                    in = inpolygon(centroid_temp_0(1),centroid_temp_0(2),contour_data{j,1}(:,1),contour_data{j,1}(:,2));
                    % check back if the centroid of the contour found belongs to
                    % the contour for the original centroid
                    if in == 1
                        centroid_temp_1 = contour_data{j,2}(1:3); % centroid of the found contour
                        % project centroid_temp_1 along y direction to xz plane
                        % where the contour_temp_0 stays.
                        in2 = inpolygon(centroid_temp_1(1),centroid_temp_1(3),contour_temp_0(:,1),contour_temp_0(:,3));
                    else
                        in2 = 0;
                    end
                end
            case 3 % contour_temp_0 is in the z direction
                if contour_data{j,2}(4) == 1 % project along x direction to yz plane
                    in = inpolygon(centroid_temp_0(2),centroid_temp_0(3),contour_data{j,1}(:,2),contour_data{j,1}(:,3));
                    % check back if the centroid of the contour found belongs to
                    % the contour for the original centroid
                    if in == 1
                        centroid_temp_1 = contour_data{j,2}(1:3); % centroid of the found contour
                        % project centroid_temp_1 along z direction to xy plane
                        % where the contour_temp_0 stays.
                        in2 = inpolygon(centroid_temp_1(1),centroid_temp_1(2),contour_temp_0(:,1),contour_temp_0(:,2));
                    else
                        in2 = 0;
                    end
                elseif contour_data{j,2}(4) == 2 % project along y direction to xz plane
                    in = inpolygon(centroid_temp_0(1),centroid_temp_0(3),contour_data{j,1}(:,1),contour_data{j,1}(:,3));
                    % check back if the centroid of the contour found belongs to
                    % the contour for the original centroid
                    if in == 1
                        centroid_temp_1 = contour_data{j,2}(1:3); % centroid of the found contour
                        % project centroid_temp_1 along z direction to xy plane
                        % where the contour_temp_0 stays.
                        in2 = inpolygon(centroid_temp_1(1),centroid_temp_1(2),contour_temp_0(:,1),contour_temp_0(:,2));
                    else
                        in2 = 0;
                    end
                elseif contour_data{j,2}(4) == 3 % project along z direction to xy plane
%                     in = inpolygon(centroid_temp_0(1),centroid_temp_0(2),contour_data{j,1}(:,1),contour_data{j,1}(:,2));
%                     % check back if the centroid of the contour found belongs to
%                     % the contour for the original centroid
%                     if in == 1
%                         centroid_temp_1 = contour_data{j,2}(1:3); % centroid of the found contour
%                         % project centroid_temp_1 along z direction to xy plane
%                         % where the contour_temp_0 stays.
%                         in2 = inpolygon(centroid_temp_1(1),centroid_temp_1(2),contour_temp_0(:,1),contour_temp_0(:,2));
%                     else
                        in2 = 0;
%                     end
                end
        end
        if in2 == 1
            contour_belong(j,:) = contour_data(j,:);
            contour_matrix(i,j) = 1;
        end
    end
    contourOnly_belong(i,:) = contour_belong(:,1); % only save contour data as a matrix of cells
    cell_inf(i,:,:) = contour_belong; % save contour data and all other informations as a two D matrix
    clear contour_belong;
end
% contour_matrix = ~isempty(contourOnly_belong);
% plot3(centroid_temp_0(1),centroid_temp_0(2),centroid_temp_0(3),'ks')
% plot3(contour_temp_0(:,1),contour_temp_0(:,2),contour_temp_0(:,3),'k-')
% plot3(contour_data{6,1}(:,1),contour_data{6,1}(:,2),contour_data{6,1}(:,3),'k-')
% plot3(contour_data{6,2}(1),contour_data{6,2}(2),contour_data{6,2}(3),'ks')

% save('Contour_all_inf.mat','contourOnly_belong','cell_inf','contour_matrix','contour_data');
hold off;

%% assign contour to cells
% add the contour itself
contour_matrix(1:1+size(contour_matrix,1):end) = 1;
sum_log = sum(contour_matrix);
bad_contour = 0;
for t = 1:length(sum_log)
    if sum_log(t) == 1 % means this contour doesn't interact with other contours, toss it out
        contour_matrix(t,t) = 0;
        bad_contour = bad_contour + 1;
%         t %(4,11,12,29 are bad contours)
    end
end

%% Using adj2cluster.m
cell_assigned_ac = adj2cluster(contour_matrix);
cells_final2 = cell(1,length(cell_assigned_ac));


figure;hold on;
for i = 1:length(cell_assigned_ac)
    cells_final2{i} = contour_data(cell_assigned_ac{i});
    cells_final_temp = cat(1,cells_final2{i});
    color_ind = rand(1,3);
    if length(cells_final_temp) > 1
        for m = 1:length(cells_final_temp)
            plot3(cells_final_temp{m}(:,1),cells_final_temp{m}(:,2),cells_final_temp{m}(:,3),'-*','color',color_ind);
        end
    else
        for m = 1:length(cells_final_temp)
            plot3(cells_final_temp{m}(:,1),cells_final_temp{m}(:,2),cells_final_temp{m}(:,3),'r-s');
        end
    end
end
title('All Cell Contours');

cell_1contour = 0;
figure;hold on;
for i = 1:length(cell_assigned_ac) 
    cells_final2{i} = contour_data(cell_assigned_ac{i});
    cells_final_temp = cat(1,cells_final2{i});
    color_ind = rand(1,3);
    if length(cells_final_temp) > 1
        for m = 1:length(cells_final_temp)
            plot3(cells_final_temp{m}(:,1),cells_final_temp{m}(:,2),cells_final_temp{m}(:,3),'-*','color',color_ind);
        end
    else
%         for m = 1:length(cells_final_temp)
%             plot3(cells_final_temp{m}(:,1),cells_final_temp{m}(:,2),cells_final_temp{m}(:,3),'r-s');
%         end
cell_1contour = cell_1contour + 1;
    end
%     hold off;
end
title('Cells with more than 1 Contour');

% save('final_assigned_cells.mat','cells_final2');

% Get the centroid information
cells_contours = cell(1,length(cell_assigned_ac));
cells_centroid = cell(1,length(cell_assigned_ac));
cells_contours_good = cell(1,length(cell_assigned_ac));
cells_centroid_good = cell(1,length(cell_assigned_ac));

for i = 1:length(cell_assigned_ac)
    cells_contours{i} = contour_data(cell_assigned_ac{i},1);
    cells_centroid{i} = contour_data(cell_assigned_ac{i},2);
    if length(cell_assigned_ac{i}) > 6 % assume there are at least 2 contours in each direction
        cells_contours_good{i} = contour_data(cell_assigned_ac{i},1);
        cells_centroid_good{i} = contour_data(cell_assigned_ac{i},2);
    else
        cells_contours_good{i} = [];
        cells_centroid_good{i} = [];
    end
end
cells_contours_good_nonEmpty = cells_contours_good(~cellfun('isempty',cells_contours_good)); % remove empty elements
cells_centroid_good_nonEmpty = cells_centroid_good(~cellfun('isempty',cells_centroid_good)); % remove empty elements

% plot cells and centroid
figure;hold on;
title('Final processed contours');
for i = 1:length(cell_assigned_ac) 
    cells_final2{i} = contour_data(cell_assigned_ac{i});
    cells_final_temp = cat(1,cells_final2{i});
    cells_centroid_temp = contour_data(cell_assigned_ac{i},2);
    color_ind = rand(1,3);
    if length(cells_final_temp) > 1
        for m = 1:length(cells_final_temp)
            plot3(cells_final_temp{m}(:,1),cells_final_temp{m}(:,2),cells_final_temp{m}(:,3),'-','color',color_ind);
            plot3(cells_centroid_temp{m}(1),cells_centroid_temp{m}(2),cells_centroid_temp{m}(3),'o','color',color_ind);
        end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
    else
%         for m = 1:length(cells_final_temp)
%             plot3(cells_final_temp{m}(:,1),cells_final_temp{m}(:,2),cells_final_temp{m}(:,3),'r-s');
%         end
% cell_1contour = cell_1contour + 1;
    end
%     hold off;
end

figure;hold on;
color_ind = rand(1,3);
for m = 1:length(cluster_1_contour)
plot3(cluster_1_contour{m}(:,1),cluster_1_contour{m}(:,2),cluster_1_contour{m}(:,3),'-','color',color_ind);
plot3(cluster_1_centroid{m}(1),cluster_1_centroid{m}(2),cluster_1_centroid{m}(3),'o','color',color_ind);
end
figure;hold on;
color_ind = rand(1,3);
for m = 1:length(cluster_1_contour)
plot3(cluster_1_centroid{m}(1),cluster_1_centroid{m}(2),cluster_1_centroid{m}(3),'o','color',color_ind);
end

cell9 = sort(cell9');
for m = 1:length(cell9)
plot3(cluster_1_centroid{cell9(m)}(1),cluster_1_centroid{cell9(m)}(2),cluster_1_centroid{cell9(m)}(3),'b*');
end
 manual_centroid{9} = cell9;

 figure;hold on;
 for m = 1:length(manual_centroid)
     cell_ind = manual_centroid{m};
     for n = 1:length(cell_ind)
%          plot3(cluster_1_contour{cell_ind(n)}(:,1),cluster_1_contour{cell_ind(n)}(:,2),cluster_1_contour{cell_ind(n)}(:,3),'r-');
         plot3(cluster_1_centroid{cell_ind(n)}(1),cluster_1_centroid{cell_ind(n)}(2),cluster_1_centroid{cell_ind(n)}(3),'r*');
     end
 end
 
figure;hold on;
for i = 1:length(cells_centroid_good_nonEmpty)
color_ind = rand(1,3);
centroid_temp = cat(1,cells_centroid_good_nonEmpty{i});
for j = 1:length(centroid_temp)
plot3(centroid_temp{j}(1),centroid_temp{j}(2),centroid_temp{j}(3),'o','color',color_ind);
end
end

% Crop data
contour_cropped = cell(1,length(cells_centroid_good_nonEmpty));
centroid_cropped = cell(1,length(cells_centroid_good_nonEmpty));
for i = 1:length(cells_centroid_good_nonEmpty)
    contour_temp = cells_contours_good_nonEmpty{i};
    centroid_temp = cells_centroid_good_nonEmpty{i};
    
    for j = 1:length(centroid_temp)
        if centroid_temp{j}(1) > 0 && centroid_temp{j}(1) < 40 && centroid_temp{j}(2) > 0 && centroid_temp{j}(2) < 40
            contour_cropped{i}{j,1} = contour_temp{j,1};
            centroid_cropped{i}{j,1} = centroid_temp{j,1}(1:3);
        end
    end
    if ~isempty(centroid_cropped{i})
    contour_cropped{i} = contour_cropped{i}(~cellfun('isempty',contour_cropped{i}));
    centroid_cropped{i} = centroid_cropped{i}(~cellfun('isempty',centroid_cropped{i}));
    end
end
contour_cropped = contour_cropped(~cellfun('isempty',contour_cropped));
centroid_cropped = centroid_cropped(~cellfun('isempty',centroid_cropped));

figure;hold on;
for i = 1:length(centroid_cropped)
color_ind = rand(1,3);
cropped_centroid_temp = cat(1,centroid_cropped{i});
cropped_contour_temp = cat(1,contour_cropped{i});
for j = 1:length(cropped_centroid_temp)
plot3(cropped_centroid_temp{j}(1),cropped_centroid_temp{j}(2),cropped_centroid_temp{j}(3),'o','color',color_ind);
plot3(cropped_contour_temp{j}(:,1),cropped_contour_temp{j}(:,2),cropped_contour_temp{j}(:,3),'-*','color',color_ind);
end
end

figure;hold on;
for i = 1:length(centroid_manual_cells)
    color_ind = rand(1,3);
    plot3(centroid_manual_cells{i}(:,1),centroid_manual_cells{i}(:,2),centroid_manual_cells{i}(:,3),'o','color',color_ind)
end

centroid_manual_cells = cell(1,18);

j = 18;
centroid_temp_3 = cursor_18;
centroid_temp_4 = nan(length(centroid_temp_3),3);

for i = 1:length(centroid_temp_3)
    centroid_temp_4(i,:) = centroid_temp_3(i).Position;
% plot3(cell1(i).Position(1),cell1(i).Position(2),cell1(i).Position(3),'o','color',color_ind)
end
centroid_manual_cells{j} = centroid_temp_4;


% %% sum interactive contours
% cell_assigned = nan(length(contour_matrix));
% for i = 1:length(contour_matrix)
%     cell_row_temp = contour_matrix(i,:);
%     cell_col_temp2 = contour_matrix(i,:);
%     for j = 1:length(cell_row_temp)
% %         if cell_row_temp(j) ~= 0 && ~isnan(cell_row_temp(j))
%         if cell_row_temp(j) ~= 0
%             cell_col_temp1 = contour_matrix(:,j);
%             cell_col_temp2 = cell_col_temp2 + cell_col_temp1';
%             for k = 1:length(cell_col_temp1)
%                 if cell_col_temp1(k) ~= 0
%                     cell_col_temp2 = cell_col_temp2 + contour_matrix(k,:);
%                 end
%             end
%         end
%     end
%     cell_col_temp2(cell_col_temp2 >= 1) = 1;
%     cell_assigned(i,:) = cell_col_temp2;
%     clear cell_col_temp2;
% end
% %% find same cells
% cell_assigned_unique = unique(cell_assigned,'row');
% cell_assigned_unique(sum(cell_assigned_unique,2) == 0,:) = [];
% cell_assigned_unique_temp = cell_assigned_unique;
% 
% cell_assigned_unique_temp(cell_assigned_unique_temp==0) = nan;
% cell_assigned_final = nan(size(cell_assigned_unique));
% f = 1;
% while ~isempty(cell_assigned_unique_temp)
%     cell_temp1 = cell_assigned_unique_temp(1,:); % assign the first row as cell_1    
%     %     num_con_cell_temp1 = sum(~isnan(cell_temp1)); % sum nonnan in cell_1
%     num_con_cell_temp1 = sum(cell_temp1); % how many contours in cell_1
% %     num_con_cell_temp2 = sum(cell_assigned_unique_temp,2);    
%     %     cell_temp1(isnan(cell_temp1)) = 0; % assign nan back to 0
%     %     cell_assigned_unique_temp(isnan(cell_assigned_unique_temp)) = 0; % assign nan back to 0    
%     diff = abs(bsxfun(@minus, cell_assigned_unique_temp, cell_temp1));
%     diff(isnan(diff)) = 1; % assign nan to 1
%     num_con_same = sum(diff == 0,2); % how many contorus are the same, means how many 0
%     cell_temp2 = cell_assigned_unique_temp(num_con_same >= 1,:); % select cells with at lest 1 same contour
%     cell_temp2(isnan(cell_temp2)) = 0; % assign nan to 0
%     cell_temp2 = sum(cell_temp2);
%     cell_temp2(cell_temp2 >= 1) = 1;
%     cell_assigned_final(f,:) = cell_temp2;
%     
%     cell_assigned_unique_temp(num_con_same >= 1,:) = [];
%     %     num_con_diff_ratio = num_con_diff./num_con_cell_temp2;
%     f = f + 1;
% end
% 
% cells_final = cell(1,length(cell_assigned_unique(:,1)));
% 
% for h = 1:length(cells_final)
%     cells_final{h} = contour_data(cell_assigned_unique(h,:) == 1,:);
% end
% 
% figure;hold on;
% for n = 1:length(cells_final)
%     cells_final_temp = cells_final{1,n};
%     color_ind = rand(1,3);
%     for m = 1:length(cells_final_temp)
%         plot3(cells_final_temp{m,1}(:,1),cells_final_temp{m,1}(:,2),cells_final_temp{m,1}(:,3),'-*','color',color_ind);
%     end
% end
% hold off;

% cell_assigned_ac2 = cell_assigned_ac(cellfun('length',cell_assigned_ac)>1);
% cell_assigned_temp = [cell_assigned zeros(length(cell_assigned),1)];
% for i = 1:length(cell_assigned)
%     cell_temp1 = cell_assigned(i,:);
%     for j = 1:length(cell_assigned)
%         cell_temp2 = cell_assigned(j,:);
%         %         if sum(abs(cell_temp1 - cell_temp2))/(length(cell_assigned) - bad_contour) < 10
%         %save these contour and then toss out these contour from
%         %'cell_assigned', so need a while loop to do this
%         if cell_temp1 == cell_temp2
%             cell_assigned_temp(i,j) = 
%         end
%     end
% end

% figure;hold on;
% for ii = 1:length(contour_data)
% plot3(contour_data{ii,1}(:,1),contour_data{ii,1}(:,2),contour_data{ii,1}(:,3),'r-*')
% end
% 
% ii=4;plot3(contour_data{ii,1}(:,1),contour_data{ii,1}(:,2),contour_data{ii,1}(:,3),'g-o')
% ii=11;plot3(contour_data{ii,1}(:,1),contour_data{ii,1}(:,2),contour_data{ii,1}(:,3),'g-o')
% ii=8;plot3(contour_data{ii,1}(:,1),contour_data{ii,1}(:,2),contour_data{ii,1}(:,3),'go')
% ii=9;plot3(contour_data{ii,1}(:,1),contour_data{ii,1}(:,2),contour_data{ii,1}(:,3),'bo')
% ii=10;plot3(contour_data{ii,1}(:,1),contour_data{ii,1}(:,2),contour_data{ii,1}(:,3),'yo')
% ii=2;plot3(contour_data{ii,2}(1),contour_data{ii,2}(2),contour_data{ii,2}(3),'ro')
% ii=8;plot3(contour_data{ii,2}(1),contour_data{ii,2}(2),contour_data{ii,2}(3),'go')
% ii=9;plot3(contour_data{ii,2}(1),contour_data{ii,2}(2),contour_data{ii,2}(3),'bo')
% ii=10;plot3(contour_data{ii,2}(1),contour_data{ii,2}(2),contour_data{ii,2}(3),'yo')


