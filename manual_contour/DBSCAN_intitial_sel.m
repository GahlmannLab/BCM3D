function [IDX, isnoise, answer1, ind_good_cluster] = DBSCAN_intitial_sel(X,epsilon,MinPts)
anser1 = 'No';
while strcmp(anser1,'No')
    % Input epsilon and MinPts for DBSCAN
    prompt = {'Epsilon','MinPts'};
    title1 = 'Parameters for DBSCAN';
    dims = [1 25];
    definput = {num2str(epsilon),num2str(MinPts)};
    in_para = inputdlg(prompt,title1,dims,definput);
    epsilon = str2double(in_para(1));
    MinPts = str2double(in_para(2));
    
    % Run DBSCAN Clustering Algorithm
    [IDX, isnoise] = DBSCAN(X,epsilon,MinPts);
    f3 = figure;
    PlotClusteringResult(X, IDX);
    
    answer1 = questdlg('Can you find a good cluster?','Yes','No');
    
    if strcmp(answer1,'Yes')
        prompt = {'Select the index of the good cluster'};
        title1 = 'Select a cluster';
        dims = [1 25];
        definput = {'1'};
        in_para = inputdlg(prompt,title1,dims,definput);
        ind_good_cluster = str2double(in_para(1));
    else
        close(f3);
    end
end
end