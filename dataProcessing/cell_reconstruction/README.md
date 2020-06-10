Run batch_process.m to reconstruct cells by building convex hull and then applying this the convex hull as mask to select volume for cells from the CNN results.

Input: Lcuts results

Output:
post_post_segments, contains voxels for each single cell.
post_post_seg_mat, save the whole processed data as a 3D matrix.

This step will call the function to check if points are including by a convex hull. The code is referred from the following link.
John D'Errico (2020). Inhull (https://www.mathworks.com/matlabcentral/fileexchange/10226-inhull), MATLAB Central File Exchange. Retrieved May 22, 2020.
