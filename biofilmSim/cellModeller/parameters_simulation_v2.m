% See detailed simulation in cellModeller https://haselofflab.github.io/CellModeller/
%
% take two .csv files 'cell_parameters.csv', 'orientation.csv'
% cell_orientation gives the angles to w,y,z axis of each cell 
% first row of the cell_parameter is cell ID, second is cell length l, 
% r is coded to be 0.4, thrid, forth and fifth are the x,y,z position of
% the center of mass
% the unit is in nanometer now
%07-17-2019 Yibo
clc;clear all;close all;

%save the parameters(positions, l and so on)
[fileName,dirName] = uigetfile('*','Load the cell parameter data',...
 'MultiSelect', 'on');
addpath(genpath(dirName));
number_of_file_parameter = length(fileName);
cell_orientation = {};
cell_parameter = {};
for i = 1: number_of_file_parameter
    cell_parameter{1,i} = readmatrix(fileName{1,i});
    for j = 1: length(cell_parameter{1,i})
        %radius = 700 + 150 * rand() ; % radomized radius in nanometers
        %cell_parameter{1,i}(2,j) = radius;
    end
end

%%
%save the orientation 
[fileName,dirName] = uigetfile('*','Load the orientation data',...
 'MultiSelect', 'on');
number_of_file_orientation = length(fileName);

for i = 1: number_of_file_orientation
    cell_orientation{1,i} = readmatrix(fileName{1,i});
    
end

%%
% concatenate to make a new file 
% each cell array is a frame
% cellid, radius, cell length l, position x,y,z, oreintation/direction(unit vector) psi x, psi y, psi z. 
parameter_eachcell = {};
for i = 1: number_of_file_parameter
    parameter_eachcell{1,i} = [cell_parameter{1,i}'  cell_orientation{1,i}];
end
save('cell_parameter_eachframe.mat','parameter_eachcell');






