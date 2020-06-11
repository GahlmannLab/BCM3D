# BCM3D
BCM3D was a image analysis workflow that combines deep learning with mathematical image analysis to accurately segment and classify single bacterial cells within biofilms in 3D fluorescence images. 

Here is a brief introduction of each module (examples and details can be found in each module separately):

1 biofilmSim: Generate 3D fluorescence images and ground truth from CellModeller results.

 Input - cell arrangement generate by CellModeller (.csv file includes cell length, radius, position, orientation).

 Output - 3D fluorescence images and corresponding ground truth.



2 trainingDataGenerate: Generate training data for CNNs. 

Input - selected ground truth

Output - labeled data for CNNs training purpose (i.e. 3 classes labeled data: 'backgound', 'cellinterior','cellboundary')



3 dataProcessing: Process output data from CNNs.

Input - confidence map generate by CNNs

Output - segment result by thresholding and segment result by further LCuts (details of LCuts can be found here: https://github.com/jwang-c/Postprocessing-using-LCuts)



4 evaluation: Evaluate segment result.

Input - segment results and ground truth (or manual annotations)

Output - evaluation result (i.e. cell counting accuracy & Jaccard index)



5 manual annotation: Code for generating manual annotation for 3D experimental data



6 utilFiles: Commonly used function. 



7 experimental data:  All experimental data that we shown in the manuscript.








