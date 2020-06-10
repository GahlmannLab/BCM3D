classdef CellModeller_Convolution_final
    properties
        name
        Distance_ratio
        SBR
        PSF
        stain
        n_of_replicate
        radius
        cellLength
        position
        orientation
        filePath
        CellModeller_data
    end
    methods
        function obj=CellModeller_Convolution_final(name,Distance_ratio,SBR,PSF,stain,...
                parameter_eachcell,n_of_replicate,filePath)
               % class constructor
                 obj.name = name;
                 obj.Distance_ratio   = Distance_ratio;
                 obj.SBR    = SBR;
                 obj.PSF  = PSF;
                 obj.stain   = stain;
                 obj.n_of_replicate = n_of_replicate;
                 obj.radius = parameter_eachcell{1,n_of_replicate}(:,2);
                 obj.cellLength = parameter_eachcell{1,n_of_replicate}(:,3);
                 obj.position = parameter_eachcell{1,n_of_replicate}(:,4:6);
                 obj.orientation = parameter_eachcell{1,n_of_replicate}(:,7:9);
                 obj.filePath = filePath;
                 obj.CellModeller_data = parameter_eachcell{1,n_of_replicate};
        end
        % get whole cell length
        function totalL = get_total_length(obj)
            totalL = 2* obj.radius + obj.cellLength ;
        end

        % show cellmodeler result and convolve it with prefered psf
        function Cellmodeller_convolutionV4 (obj, psf, n)
            %% Load parameter and catagorize
            % loading experimental psf
            [psf_conv, psf_decon, background] = LoadPSF_select(psf);
            %creating biofilm volume
            %CellModeller_data = parameter_eachcell{1,obj.n_of_replicate};
            rng(10);
            simbox=[23000 23000 13000]; %whole volume % the Z has to be bigger than 101 (the z of PSF)
            simbox=ceil(simbox.*obj.Distance_ratio);
            voxelSize=100;%voxel size
            [numCells,rawLocs,ground,surface_select_points,...
                whole_select_points,ground_label]=cellVolume_cellmodeller(obj.CellModeller_data,...
                voxelSize,simbox,obj.Distance_ratio,n);
            % output membrane and cytosolic points
            if strcmp(obj.stain,'cytosol') == 1
                select_points = whole_select_points;

            elseif strcmp(obj.stain,'membrane') == 1
                select_points = surface_select_points;
            else
                error('no matched stain');
            end
            str1 = num2str(numCells);
            str7 = num2str(obj.Distance_ratio);
            number = num2str(obj.n_of_replicate);

            %% convolution
            model_data = convn(select_points,psf_conv,'same');%convolution simulate model with psf

            % simulate noise and background
            vari = 3.04;% Gaussian-distributed camera read-out noise
            noise = poissrnd(background,size(model_data)) + normrnd(zeros(size(model_data)),vari);
            % Sum possion noise and gaussina noise
            %% add noise and set SBR
            str8 = num2str(obj.SBR);
            SignalPerCell=sum(model_data(:))/numCells;
            SBR_initial = SignalPerCell/background; %photons per cell/phonon per voxel
            model_data_temp = model_data.*(( obj.SBR-1)/(SBR_initial-1)); % change signal level to setting SBR
            output_data = model_data_temp + noise;
            %% deconvolve
            rawdata = single(output_data)-background;%%
            rawdata(rawdata<0) = 0;% ensure signal is larger than 0
            nIter=10;
            rawdata_size=size(rawdata);
            psf_decon_rot_size=size(psf_decon);

            if psf_decon_rot_size > rawdata_size
                %need pad in z direction for further deconv
                rawdata(:,:,rawdata_size(3)+1:psf_decon_rot_size(3))=0;
            end

            deconvolved = deconvlucy(rawdata, psf_decon, nIter) * numel(rawdata);% decovolution here, after deconv, intensity is very low (not sure why), so multiply a large number to enhance intensity
            deconvolved = deconvolved(:,:,1:rawdata_size(3));
            %% output data name
            %str2 = '_cells_without_noise_background';
            str3 = '_cells_with_noise_background';
            str4 = '_deconvolved';
            %str5 = '_cells_ground_truth';
            %str6 = '_cells_ground_truth_annotation.tif';
            %% specify or make output directories 
            formatOut = 'mmddyyyy';
            folderDate = datestr(now,formatOut);
            decon_image_dir = strcat('./result/',obj.stain,'/',folderDate,'/','deconvolved_image','/');
            label_dir = strcat('./result/',obj.stain,'/',folderDate,'/','ground_truth','/');
            raw_dir = strcat('./result/',obj.stain,'/',folderDate,'/','raw_img','/');
            if ~exist(decon_image_dir, 'dir')
               mkdir(decon_image_dir);
            end
            if ~exist(label_dir, 'dir')
               mkdir(label_dir);
            end
            if ~exist(raw_dir, 'dir')
               mkdir(raw_dir);
            end
            raw_image_name = strcat(raw_dir,number,'_SBR_',str8,'_DR_', ...
                str7, '_', obj.stain, '_', str1, str3);
            label_image_name = strcat(label_dir,number,'_','SBR_',str8,'_DR_', str7,'_',...
                str1, str4, obj.stain, '_Label');
            decon_image_name=strcat(decon_image_dir, number,'_','SBR_',str8,'_DR_', str7,'_',...
                str1, str4, obj.stain, '_T1');
            %% write .nii files for Niftynet and tif files 
            % output raw images
            write3Dtiff_V2(uint16(output_data),strcat(raw_image_name,'.tif'));
                      
            % output deconvolved .nii images
            seg_img = imrotate(uint16(deconvolved), 90);
            seg_img = flip(seg_img, 1);
            niftiwrite(seg_img, strcat(decon_image_name, '.nii'));

            % output .tif label files
            write3Dtiff_V2(uint16(ground_label),strcat(label_image_name, '.tif'));
        end
    end
end



