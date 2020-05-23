1. Run batch_process.m to reconstruct cells by building convex hull and then applying this the convex hull as mask to select volumn for cells from the CNN results.
This step will call the function to check if points are including by a convex hull. The code is referred from the following link.
John D'Errico (2020). Inhull (https://www.mathworks.com/matlabcentral/fileexchange/10226-inhull), MATLAB Central File Exchange. Retrieved May 22, 2020.
