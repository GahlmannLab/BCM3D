#### Calculate evaluation metrics

1. calculate Jacard Index and cell counting accuracy (also called Average Precision)

   plot_evaluation_metrics.m also calculates Dice, False positive percentage,
   False negative percentage, Distance of major axis, and relative angle between
   major axis. 
   When each cell  matches with  the GT cell, it  gets deleted from the list. There
   can't be cells that are matched twice.
    *inputs*: have to be a  directory that contains tif ground truth images and a directory that contains nii prediction masks from Unet.

2. calculate local density based on tiling of the image.

   calc_local_density_v3_tiles.m

   Tile GT images/3D mannual annotations and calculate local density:
   local density = cell volume / tile volume.
   *two output metrics*: maximum local density and the mean of top 10 local densities.

3. calculate the mean and median  intercellular distances.

   calc_intercellular_distances.m
   
   *input*: GT or segmentation results.
   step 1: get all pixels that are on the boundary of each object.
   step 2: iteratively search the minimum distance for each pixel within an object using knnsearch. Get a minimum distance for that object. 
   step 3: Go over all objects to find minimum distances distribution for all objects.
   *outputs*: meanD_results, medianD_results, and minDistTable for all distances.

