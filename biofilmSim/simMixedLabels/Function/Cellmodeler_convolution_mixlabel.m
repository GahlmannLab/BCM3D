%% Generate 3D fluorescence stack of mixed labeled biofilms
% (membrane and membrane&cytosolic labeled) based on CellModeller simulation 


%% Loading experimental PSF (point spread function)
[psf_rot, psf_decon_rot, background_1] = LoadPSF();

%% Creating biofilm volume
rng(10);
simbox=[22000 22000 12000]; %Size of volume
[numCells,surf_num, surf_interior_num, surface_select_points,surf_interior_select_points,ground_surf_label, ground_surf_interior_label]=cellVolume_cellmodeler_mixlabel(CellModeller_data,voxelSize,simbox,s_percent);

%% Preallocate data
mix_data = zeros(size(surface_select_points));
select_points = {surface_select_points, surf_interior_select_points};
ground_label = {ground_surf_label, ground_surf_interior_label};
%% Simulate noise and background
vari = 3.04; % Gaussian-distributed camera read-out noise
noise = poissrnd(background_1,size(mix_data)) + normrnd(zeros(size(mix_data)),vari); 
% Sum possion noise and gaussina noise

%% Convolve with PSF
type_name = {'surf', 'surfinterior'};
for n =1:length(select_points)
model_data = convn(select_points{n},psf_rot,'same');

%% Calculate initial SBR
if n==1 && surf_num>0%Membrane labeled cells
SignalPerCell=sum(model_data(:))/surf_num;
SBR_initial = SignalPerCell/background_1; %signal per cell/background
model_data_temp = model_data.*((SBR-1)/(SBR_initial-1)); 
% change signal level to setting SBR  
datafolder = datafolder3;%datafolder3 is for membrane labeled cells
elseif n==2 && surf_interior_num > 0
% membrane&cytosolic labeled cells, SBR should be twice larger
SignalPerCell=sum(model_data(:))/surf_interior_num;
SBR_initial = SignalPerCell/background_1;
SBR = 2*SBR;
model_data_temp = model_data.*((SBR-1)/(SBR_initial-1)); % Change signal level to setting SBR
datafolder =datafolder4;%datafolder4 is for membrane&cytosolic labeled cells
end

%% Get mix data
mix_data = mix_data + model_data_temp;

%% Output data
str1 = num2str(numCells);

filename5=strcat(datafolder,prefix, str1, '_',type_name{n},'_Label.tif');%

write3Dtiff_V2(uint16(ground_label{n}),filename5)% Output and save groundtruth
end
%% Add noise to mix data
mix_data = mix_data+noise;
write3Dtiff_V2(uint16(mix_data),strcat(datafolder1,prefix,str1 ,'_raw.tif')) 
% Output and save raw data
%% Deconv mix data
rawdata = single(mix_data)-background_1;
rawdata(rawdata<0)=0;
nIter=10;
deconvolved = deconvlucy(rawdata, psf_decon_rot, nIter) * numel(rawdata);
niftiwrite_v2(uint16(deconvolved),strcat(datafolder2,prefix,str1 ,'_deconv_T1.nii')); 
% Output and save deconv data