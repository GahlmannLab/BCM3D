function [maxLocalDensity, mean10LocalDensity] = calc_tile_maximum_density(img,tile,incrementxy,incrementz )
%% devide image into local cubes 
    %parameters
%     tile = [64 64 8];
%     incrementxy = 8; %xy overlaps
%     incrementz = 4; % z overlaps

    kernalVolume = tile(1) * tile(2) * tile(3);
    blockSizeRow = tile(1);blockSizeCol = tile(2);blockSizeZ = tile(3);

    %% cal density of each block
    % blockSizeRow: Rows in block
    % blockSizeCol: Columns in block
    [nrows, ncols, nz] = size(img);

    % Calculate size of each block by rows and columns
    n = 1;
    for i=1:incrementxy:nrows-blockSizeRow+1
         for j = 1:incrementxy:ncols-blockSizeCol+1
             for k = 1:incrementz:nz-blockSizeZ+1
                tileImages = img(i:i+blockSizeRow-1,...
                    j:j+blockSizeCol-1,...
                    k:k+blockSizeZ-1); 

                 if sum(tileImages,'all')== 0
                    continue
                    else   
                        cellVolume = get_cellVolume(tileImages);
                        density(n,1) =  cellVolume/kernalVolume * 100;
                 end
                 n = n +1;
             end
          end
    end

    %density_noOutlier = rmoutliers(density(density~= 0));
    maxLocalDensity = max(density);
    mean5LocalDensity = mean(maxk(density,5));
    mean10LocalDensity = mean(maxk(density,10));
end

function cellVolume = get_cellVolume(GT_img)
    sz = size(GT_img);
    l_gt = unique(GT_img);
    l_gt(l_gt == 0)=[];
    gt_index = cell(1,length(l_gt));
    cellVolume = 0;
    for i=1:length(l_gt)
        % calculate density of each cell and add to get total cell volume
        r1 = find(GT_img == l_gt(i));
        gt_index{i} = r1; 
        cellVolume = cellVolume+ length(r1);
    end
end
