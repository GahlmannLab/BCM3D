% manually select cell contour 
% Mingxing Zhang, Gahlmann lab, Chemistry department, University of Virginia

%% Load deskewed or deconvolved data (3D data set)
[fileName,dirName] = uigetfile('*','Load the deskewed or deconvolved data');
dataFile = [dirName fileName];
dataFileInfo = imfinfo(dataFile);
numFrames = length(dataFileInfo);
% frame_gap = 3; % select data every 'frame_gap' frames
processed_data = zeros(dataFileInfo(1).Height,dataFileInfo(1).Width,numFrames);

folderName = fileName(1:end - 4);
mkdir(folderName);
cd (folderName);

for frame = 1:numFrames
    processed_data(:,:,frame) = double(imread(dataFile,frame,'Info',dataFileInfo));
end

% select data every 3 frames
% processed_data = processed_data(1:5:dataFileInfo(1).Height,1:5:dataFileInfo(1).Width,1:3:numFrames);

% resize the 3D data set
% qx = 20; qy = 20; qz = 10;
% [x y z] = ndgrid(linspace(1,size(processed_data,1),qx),...
%     linspace(1,size(processed_data,2),qy),linspace(1,size(processed_data,3),qz));
% processed_data = interp3(processed_data,x,y,z);

scrsz = get(0,'ScreenSize');
%% Select cell contour from x, y and z direction
[x_dim, y_dim, z_dim] = size(processed_data); % the x and y dimension is reversed when import the data

% how many frames you want to select in each dimension?
prompt = {'x skip','y skip','z skip'; 'x0','y0','z0'; 'x end','y end','z end'};
title1 = 'Parameters to select frames';
dims = [1 25];
% definput = {'20','20','10'}; % finally will select [21 21 11] frames
definput = {'1','1','1'; '1','1','1'; num2str(y_dim),num2str(x_dim),num2str(z_dim)}; % finally will select [21 21 11] frames
para_frames = inputdlg(prompt,title1,dims,definput);
para_frames = str2double(para_frames);
Num_x_frames = length(para_frames(2):para_frames(1):para_frames(3));
Num_y_frames = length(para_frames(5):para_frames(4):para_frames(6));
Num_z_frames = length(para_frames(8):para_frames(7):para_frames(9));

% select cell contour from x direction, the x direction here is the same
% the x direction of the real data set
if para_frames(1) ~= 0
    loop_ind = 1;
    for k = para_frames(2):para_frames(1):para_frames(3)
        h = figure;
        %     axis square;
%         set(h,'Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1024 768]);
        hold on;
        title(['Frame ' num2str(loop_ind) 'of ' num2str(Num_x_frames) ' frames']);
        
        x_temp = squeeze(processed_data(:,k,:)); % select frames in x axis
        imagesc(x_temp);
        num_cell = 1; % to mark how many cells has been selected
        select = 1; % indicator to determine if select more cells       
        [BW,xi,yi] = roipoly;
        
        answer1 = questdlg('Would you like save this contour?','Yes','No');
        switch answer1
            case 'Yes'
                x_frames(k).cellContour{num_cell} = [xi yi];
                plot(xi,yi,'-k');
            case 'No'
                x_frames(k).cellContour{num_cell} = [];
        end
        
        while select ~= 0
            answer2 = questdlg('Would you like select more cells?','Yes','No');
            switch answer2
                case 'Yes'
                    select = 1;
                    num_cell = num_cell + 1;
                    [BW,xi,yi] = roipoly;
                    answer3 = questdlg('Would you like save this contour?','Yes','No');
                    switch answer3
                        case 'Yes'
                            x_frames(k).cellContour{num_cell} = [xi yi];
                            plot(xi,yi,'-k');
                        case 'No'
                            x_frames(k).cellContour{num_cell} = [];
                    end
                case 'No'
                    select = 0;
            end
        end
        saveas(h,['x_frames_' num2str(k) '.fig']);
        saveas(h,['x_frames_' num2str(k) '.tif']);
        save(['x_frames_' num2str(k) '.mat'], 'x_frames');
        hold off;
        close (h);
        loop_ind = loop_ind + 1;
    end
else
end

% select cell contour from y direction, the y direction here is the same
% the y direction of the real data set
if para_frames(4) ~= 0
    loop_ind = 1;
    for k = para_frames(5):para_frames(4):para_frames(6)
        h = figure;
        %     axis square;
        set(h,'Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1024 768]);
        hold on;
        title(['Frame ' num2str(loop_ind) 'of ' num2str(Num_y_frames) ' frames']);
        
        y_temp = squeeze(processed_data(k,:,:)); % select frames in x axis
        imagesc(y_temp);
        num_cell = 1; % to mark how many cells has been selected
        select = 1; % indicator to determine if select more cells        
        [BW,xi,yi] = roipoly;
        
        answer1 = questdlg('Would you like save this contour?','Yes','No');
        switch answer1
            case 'Yes'
                y_frames(k).cellContour{num_cell} = [xi yi];
                plot(xi,yi,'-k');
            case 'No'
                y_frames(k).cellContour{num_cell} = [];
        end
        
        
        while select ~= 0
            answer2 = questdlg('Would you like select more cells?','Yes','No');
            switch answer2
                case 'Yes'
                    select = 1;
                    num_cell = num_cell + 1;
                    [BW,xi,yi] = roipoly;
                    answer3 = questdlg('Would you like save this contour?','Yes','No');
                    switch answer3
                        case 'Yes'
                            y_frames(k).cellContour{num_cell} = [xi yi];
                            plot(xi,yi,'-k');
                        case 'No'
                            y_frames(k).cellContour{num_cell} = [];
                    end
                case 'No'
                    select = 0;
            end
        end
        saveas(h,['y_frames_' num2str(k) '.fig']);
        saveas(h,['y_frames_' num2str(k) '.tif']);
        save(['y_frames_' num2str(k) '.mat'], 'y_frames');
        hold off;
        close (h);
        loop_ind = loop_ind + 1;
    end
else
end

% select cell contour z direction, the z direction here is the same
% the z direction of the real data set
if para_frames(7) ~= 0
    loop_ind = 1;
    for k = para_frames(8):para_frames(7):para_frames(9)
        h = figure;
        %     axis square;
        set(h,'Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1024 768]);
        hold on;
        title(['Frame ' num2str(loop_ind) 'of ' num2str(Num_z_frames) ' frames']);
        
        z_temp = squeeze(processed_data(:,:,k)); % select frames in x axis
        imagesc(z_temp);
        num_cell = 1; % to mark how many cells has been selected
        select = 1; % indicator to determine if select more cells
        
        [BW,xi,yi] = roipoly;
        
        answer1 = questdlg('Would you like save this contour?','Yes','No');
        switch answer1
            case 'Yes'
                z_frames(k).cellContour{num_cell} = [xi yi];
                plot(xi,yi,'-k');
            case 'No'
                z_frames(k).cellContour{num_cell} = [];
        end
        
        while select ~= 0
            answer2 = questdlg('Would you like select more cells?','Yes','No');
            switch answer2
                case 'Yes'
                    select = 1;
                    num_cell = num_cell + 1;
                    [BW,xi,yi] = roipoly;
                    answer1 = questdlg('Would you like save this contour?','Yes','No');
                    switch answer3
                        case 'Yes'
                            z_frames(k).cellContour{num_cell} = [xi yi];
                            plot(xi,yi,'-k');
                        case 'No'
                            z_frames(k).cellContour{num_cell} = [];
                    end                    
                case 'No'
                    select = 0;
            end
        end
        saveas(h,['z_frames_' num2str(k) '.fig']);
        saveas(h,['z_frames_' num2str(k) '.tif']);
        save(['z_frames_' num2str(k) '.mat'], 'z_frames');
        hold off;
        close (h);
        loop_ind = loop_ind + 1;
    end
else
end
% save('cell_contour.mat','x_frames','y_frames','z_frames')


