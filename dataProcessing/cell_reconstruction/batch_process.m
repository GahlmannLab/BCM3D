% load lcuts results
[lcuts_files, lcuts_file_path] = uigetfile({'*.mat'},...
    'Open LCuts results',...
    'MultiSelect', 'on');

for j = 1:length(lcuts_files)
    lcuts_file_path_temp = fullfile(lcuts_file_path, lcuts_files{j}); %get input_tiff with absolute path
    file_name_temp = lcuts_files{j};
    file_name_temp = file_name_temp(1:end-4);
    newFileName = strcat(file_name_temp, '_pp');
    LcutsBendingCells_bwdist(lcuts_file_path_temp, newFileName);
    display(lcuts_files{j});
end