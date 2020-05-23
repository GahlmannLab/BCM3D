1. Run manual_select_contour.m to manually trace cell outlines.
2. Run assign_contour2cells_v2.m to group 2D cell contours to single cells
This step will call the function to extract clusters from adjecent matrix. The code is referred from the following link.
Raphaël Candelier (2020). Adjacency matrix to clusters (https://www.mathworks.com/matlabcentral/fileexchange/60676-adjacency-matrix-to-clusters), MATLAB Central File Exchange. Retrieved May 22, 2020.

3. Run Clustering_1st.m to segment cells from unsegmented clusters in step 2.
This step will call the function DBSCAN.m to group points. The code is referred from the following link.
Yarpiz (2020). DBSCAN Clustering Algorithm (https://www.mathworks.com/matlabcentral/fileexchange/52905-dbscan-clustering-algorithm), MATLAB Central File Exchange. Retrieved May 22, 2020.