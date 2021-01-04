#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 21:27:04 2020

@author: lemuel
"""
from granules.simulation import dcdReader
from granules.structure.LAMMPSdata import LammpsData, AtomsFull, MolecularTopologyData, CharmmForceField, Box
from granules.structure.NAMDdata import NAMDdata
import matplotlib.pyplot as plt
import pandas as pd

def bondLenghtsTrajectory(lammpsData, trayectory,bondTypes=None):
    '''Function that selects the bondType of the trajectory and returns a 
       panda table with the bond lengths '''

    firstFrame = lammpsData.atomprop.atoms.copy()  # will be restored before method returns
    if bondTypes != None:        
        bondLgth = lammpsData.bondLength().loc[lammpsData.bondLength()["bType"].isin(bondTypes)]
    else:
        bondLgth = lammpsData.bondLength()
    for X,Y,Z,pbc in trayectory:
        lammpsData.atomprop.atoms.update_corrdinates((X,Y,Z))
        if bondTypes != None:          
            bondLgth = bondLgth.append(lammpsData.bondLength().loc[lammpsData.bondLength()["bType"].isin(bondTypes)])
        else:
            bondLgth = bondLgth.append(lammpsData.bondLength())
            
    lammpsData.atomprop.atoms = firstFrame
    return bondLgth

def describeBondTrajectory(lammpsData, trayectory,bondType=None):
    '''Function that describes the bonds of the selected trajectory in a panda dataframe'''
       
    return bondLenghtsTrajectory(lammpsData, trayectory,bondType).groupby('bType').describe()


def describeAngleTrajectory(lammpsData, trayectory):
    '''Function that describes the angles of a trajectory in a panda dataframe'''
    firstFrame = lammpsData.atomprop.atoms.copy()  # will be restored before method returns
    angleLgth = lammpsData.angleLength()
    for X,Y,Z,pbc in trayectory:
        lammpsData.atomprop.atoms.update_corrdinates((X,Y,Z))
        angleLgth = angleLgth.append(lammpsData.angleLength())

    lammpsData.atomprop.atoms = firstFrame
    return angleLgth.groupby('anType').describe()
    

def bondTypeHist(lammpsData, trajectory, bondTypes=None, **kwargs):
    '''Draws histogram for each of the bond types and returns the axes of
    the histograms and a DataFrame with the bond lengths of the selected bondtypes
    '''
    blt = bondLenghtsTrajectory(lammpsData, trajectory, bondTypes)
    result = {name:pd.Series() for name in blt.bType.unique()}
    for name,group in blt.groupby('bType'):
        result[name] = result[name].append(group.dist.reset_index(drop=True))
    result = pd.DataFrame(result)
    ax = result.hist(**kwargs)
    return ax, result


if __name__ == "__main__":
    SELECTED_BOND_TYPES = [5,49]
    SELECTED_BOND_TYPES = None
    #Prepares the data        
    ld = LammpsData()#atoms=AtomsFull(), topology=MolecularTopologyData(), forceField=CharmmForceField(), region=Box())
    namdChignolin = NAMDdata()
    namdChignolin.readFiles('../tests/chignolin/chignolin.pdb', '../tests/chignolin/chignolin.psf', '../tests/chignolin/chignolin.prm')
    
    ld.loadNAMDdata(namdChignolin)
    #bondLgth = ld.bondLength()
    
    trayectoria = dcdReader.Trajectory("../tests/chignolin/chignolin.dcd")
    print("Output de la tabla para el histograma:\n",describeBondTrajectory(ld, trayectoria,SELECTED_BOND_TYPES))

    trayectoria = dcdReader.Trajectory("../tests/chignolin/chignolin.dcd")
    ax, bTable = bondTypeHist(ld, trayectoria,bins=23,layout=(20,3),figsize=(10,60),bondTypes=SELECTED_BOND_TYPES)
    
    
    #print("Output del histograma:\n",bondLenghtsTrajectory(ld, trayectoria,[5,6]))
    #trayectoria = dcdReader.Trajectory("../tests/chignolin/chignolin.dcd")
    #print("Output del histograma:\n",bondTypeHist(ld, trayectoria,5))
    #print(bondLenghtsTrajectory(ld, trayectoria,5).groupby("bType"))
    #print("Output del histograma:\n",bondTypeHist(ld, trayectoria,5))
    #fig = plt.subplots((1,5))
    #result = bondLenghtsTrajectory(ld, trayectoria,[5,6])
    #r=bondTypeHist(ld, trayectoria,[5])
    
    
    #print("Output de la tabla para el histograma:\n",describeAngleTrajectory(ld, trayectoria))
    
    
    
    

   





