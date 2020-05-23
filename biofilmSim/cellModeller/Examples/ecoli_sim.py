import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=True, gamma = 10, max_planes=6)

    biophys.addPlane((0,0,0),(0,0,1),1.0) #Base plane
    biophys.addPlane((10,0,0),(-1,0,0),1.0)
    biophys.addPlane((-10,0,0),(1,0,0),1.0)
    biophys.addPlane((0,10,0),(0,-1,0),1.0)
    biophys.addPlane((0,-10,0),(0,1,0),1.0)
    biophys.addPlane((0,0,10),(0,0,-1),1.0)


    # use this file for reg too
    regul = ModuleRegulator(sim, sim.moduleName)
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(-5,-1,0.5), dir=(-1,0,0), rad=0.4)
    sim.addCell(cellType=0, pos=(6,8,0.5), dir=(-0.5,0.2,0), rad=0.4)
    sim.addCell(cellType=0, pos=(-1,-7,1), dir=(0,0,0), rad=0.4)
    #sim.addCell(cellType=0, pos=(1,11,0.5), dir=(1,0.4,0), rad=0.4)
    #sim.addCell(cellType=0, pos=(-1,7,0.5), dir=(1,1,0), rad=0.4)

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sim.pickleSteps = 10

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = 2 + random.uniform(0.0,1.5)
    # Specify growth rate of cells
    cell.growthRate = 0.2
    cell.color = (1.0,1.0,1.0)

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.iteritems():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 2 + random.uniform(0.0,1.5)
    d2.targetVol = 2 + random.uniform(0.0,1.5)
