#### CellModeller data processing

See detailed biofilm simulation in cellModeller https://haselofflab.github.io/CellModeller/.
load_pickle.ipynb takes the saved pickle data from cellModeller and saves both the cell positions and orientations in two separate csv files. The example pickle data is shown in <u>Examples</u>.

parameters_simulation_v2 takes two csv files 'cell_parameters.csv' and 'orientation.csv'.
cell_orientation gives the angles to x,y,z axis of each cell.

The first row of the cell_parameter is cell ID;second is cell length l; the third, forth and fifth are the x,y,z position of the center of mass.  The radius r, which is the radius of the spherical cylinder is either hard coded to be 400 nanometers or random numbers between biological reasonable range.

The unit is in nanometer.