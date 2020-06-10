%% Generate simulated 3D fluorscence stack of mixed shape biofilms 
%  based on CellModeller simulation
 

%% Creating biofilm volume
rng(10);
simbox=[22000 22000 12000]; %Size of volume
[numCells,rod_num, sphere_num, rod_select_points,sphere_select_points,ground_rod_label, ground_sphere_label]=cellVolume_cellmodeler_mixshape(CellModeller_data,voxelSize,simbox,r_percent);

mix_data = zeros(size(rod_select_points));% Preallocate data
select_points = {rod_select_points, sphere_select_points};
ground_label = {ground_rod_label, ground_sphere_label};

%% Simulate noise and background
vari = 3.03; %Gaussian-distributed camera read-out noise
noise = poissrnd(background,size(mix_data)) + normrnd(zeros(size(mix_data)),vari);
% Sum possion noise and gaussina noise

%% Convolve data with PSF
type_name = {'rod', 'sphere'};
for n =1:length(select_points)
model_data = convn(select_points{n},psf_conv,'same');%convolution simulate model with psf

%% Calculate initial SBR
if n==1 && rod_num>0 %Rod shaped cell
SignalPerCell=sum(model_data(:))/rod_num;
SBR_initial = SignalPerCell/background; %Signal per cell/background
model_data_temp = model_data.*((SBR-1)/(SBR_initial-1)); % Change signal level to setting SBR  
datafolder = datafolder3;% datafolder3 is for rod shaped cells

elseif n==2 && sphere_num > 0% Sphere shaped cell
SignalPerCell=sum(model_data(:))/sphere_num;
SBR_initial = SignalPerCell/background; %signal per cell/background 
SBR = 0.4*SBR; % Reduce signal of sphere shaped cell otherwise it will be much brighter than rod shaped cells (due to the different cell volume)
model_data_temp = model_data.*((SBR-1)/(SBR_initial-1)); % Change intensity to setting SBR
datafolder =datafolder4;%datafolder4 is for sphere shaped cells
end

%% get mix data
mix_data = mix_data + model_data_temp;

%% Output data
str1 = num2str(numCells);

filename5=strcat(datafolder,prefix, str1, '_',type_name{n},'_Label.tif');

write3Dtiff_V2(uint16(ground_label{n}),filename5) %Output and save groundtruth
end
%% Add noise to mix data
mix_data = mix_data+noise;
write3Dtiff_V2(uint16(mix_data),strcat(datafolder1,prefix,str1 ,'_raw.tif')) 
% Output and save raw data
%% Deconv mix data
rawdata = single(mix_data)-background;
rawdata(rawdata<0)=0;
nIter=10;
deconvolved = deconvlucy(rawdata, psf_decon, nIter) * numel(rawdata);
niftiwrite_v2(uint16(deconvolved),strcat(datafolder2,prefix,str1 ,'_deconv_T1.nii')); 
% Output and save deconv data