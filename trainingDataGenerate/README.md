#### To generate training data for CNNs

We used a common 3-classes (5-classes for heterogenous biofilms) labeling strategies in our CNNs training: label background voxels as zero, cell interiors as one and cell boundaries as two.

These scripts shown here can be used to generate 3D labels for training (NiftyNet platform). NiftyNet is a open source CNN platforms. For more details, please go to the official website, https://niftynet.io/Training. 

 Jupyter Notebook scripts:

package dependencies:  nibabel,numpy, scikit-image

All the inputs should be 3D .tif ground truth images that have each cell labeled as unique integers.
trainingDataGenerate generates labels for single-population biofilms(label background voxels as zero, cell interiors as one and cell boundaries as two.), whereas trainingDataGenerate_mixLabel and trainingDataGenerate_mixShape work for two-populations respectively (label background voxels as zero, cell interiors of population 1 as one, cell boundaries of population 1 as two, cell interiors of population 2 as three and cell boundaries of population 2 as four).
trainingDataGenerate_mixLabel: input data should have GT pairs that have prefixes of "surf__" and "surfInterior_".

 trainingDataGenerate_mixShape: input data should have GT pairs that have suffix of "rod_Label" and "sphere_Label".

Examples of pretrained models are stored in folder 'pretrainedModels', each model folder contains 5 subfolders, 'config_file' stores parameters to train the network and inference data by the trained network, 'trainingdata' stores training pairs to train the network,  'testdata' stores example data for testing, 'output' stores output of the network for the test data. '..._model' stores the network.