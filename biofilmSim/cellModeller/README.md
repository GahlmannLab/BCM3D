#### CellModeller data generation and processing

CellModeller is a Python-based framework for modelling large-scale multi-cellular systems that was developed by the Haseloff lab.

We have only tested CellModeller on Mac Mojave. For more details and installation instructions of CellModeller, see the official website https://haselofflab.github.io/CellModeller/. An example python script to use CellModeller is shown in <u>Examples</u>. CellModeller saves pickled biofilm cell arrangement data, and we used **load_pickle.ipynb** to save both the cell positions and orientations in two separate csv files. 

**parameters_simulation_v2** combines two csv files, 'cell_parameters.csv' and 'orientation.csv', to a .mat file, 'cell_parameter_eachframe.mat' for downstream biofilm simulations. The file contains cell radius, length, positions and orientations. (example .csv files are in ./csvFiles/ecoli_csv)

