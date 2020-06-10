#### Thresholding

**Thresholding** takes CNN confidence maps and output segmentation masks. We found a high threshold between 0.88 and 0.94 yield the best segmentation results with well-trained CNN models. Since cell boundaries are eroded in the CNN processing, we also prefer to dilate each objects  by two voxels. 

*inputs*: CNN confidence maps (.nii), threshold values (default is 0.94 for cytosol-labeled cells and 0.88 for membrane-labeled cells).

*outputs*: segmentation masks.