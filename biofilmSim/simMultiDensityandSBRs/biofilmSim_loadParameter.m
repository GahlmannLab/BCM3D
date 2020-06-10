 %load cell arrangment parameters and create object of CellModeller_Convolution_Class
 % 'cell_parameter_eachframe.mat' is the cell arrangement information for
 % each frame derived from CellModeller
 % 'parameter_config.csv' is the information for simulation, containing the
 % density, SBRs and labeling methods parameters
 % PSFs are experimentally measured LLSM PSFs either in 488 or 561
 % channels, both the raw PSFs and deconvolved PSFs are needed for
 % simulation
 % The outputs are raw fluorecence images, GTs, and deconvolved images of
 % the corresponding simulated biofilms
 %updated in 20200519
 
addpath(genpath('functions'));addpath(genpath('PSFs'));
addpath(genpath('../../utilFiles'));
% parameter_config.csv 
[file1, path1] = uigetfile('*.csv','Select a parameter_config csv file');
addpath(genpath(path1));
object_config = readtable(file1);
 
% Select an .mat file that contains cell arrangment information
[file2, path2] = uigetfile('*.mat','Select an .mat file that contains cell arrangment information');
addpath(genpath(path2));
fullFile2 = strcat(path2, file2);
load(fullFile2);
filePath = pwd;

% choose a psf to convolve and deconvolvechoose the
% right PSF. Options(string): 561psf and 488psf
% keep the names of the matfiles the same
prompt = {'Enter your PSF (561psf or 488psf)', 'number of emitters (default:(myxo:n = 2000, ecoli:n = 500))'};
dlgtitle = 'Chosen PSF and emitters';
definput = {'488psf', '500'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,[1 80],definput,opts);
psf_chosen = answer{1}; number_of_emitters = str2double(answer{2});

% convolution
parfor k = 1:height(object_config)
    tic
    a = object_config(k,:);
    obj=CellModeller_Convolution_final(cell2mat(a.Var1),a.Var2,a.Var3, ...
        cell2mat(a.Var4),cell2mat(a.Var5),parameter_eachcell,...
        a.Var7,filePath); 
    % run the main class function Cellmodeler_convolutionV4.m, 
    Cellmodeller_convolutionV4(obj, psf_chosen, number_of_emitters);
    disp(k);
    toc
 end 


