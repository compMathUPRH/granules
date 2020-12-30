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

def bondLenghtsTrajectory(lammpsData, trayectory,bondType=None):
    '''Function that selects the bondType of the trajectory and returns a panda table'''
    
    if bondType != None:        
        bondLgth = lammpsData.bondLength().loc[lammpsData.bondLength()["bType"].isin(bondType)]
    else:
        bondLgth = lammpsData.bondLength()
        
    for X,Y,Z,pbc in trayectory:
        lammpsData.atomprop.atoms.update_corrdinates((X,Y,Z))
        if bondType != None:          
            bondLgth = bondLgth.append(lammpsData.bondLength().loc[lammpsData.bondLength()["bType"].isin(bondType)])
        else:
            bondLgth = bondLgth.append(lammpsData.bondLength())
    return bondLgth

def describeBondTrajectory(lammpsData, trayectory,bondType=None):
    '''Function that describes the bonds of the selected trajectory in a panda dataframe'''
       
    return bondLenghtsTrajectory(lammpsData, trayectory,bondType).groupby('bType').describe()


def describeAngleTrajectory(lammpsData, trayectory):
    '''Function that describes the angles of a trajectory in a panda dataframe'''
    angleLgth = lammpsData.angleLength()
    for X,Y,Z,pbc in trayectory:
        lammpsData.atomprop.atoms.update_corrdinates((X,Y,Z))
        angleLgth = angleLgth.append(lammpsData.angleLength())

    return angleLgth.groupby('anType').describe()
    

def bondTypeHist(lammpsData, trajectory, bondType=None):
    '''Return a list of tuples, with the selected bondtype and the axes of 
    the histogram '''
    
    result = bondLenghtsTrajectory(lammpsData, trajectory, bondType)
    print(result)
    result.plot(kind='hist', column='dist', by='bType')
    return result

def bondTypeHist3(lammpsData, trajectory, bondType=None):
    '''Return a list of tuples, with the selected bondtype and the axes of 
    the histogram '''
    result = pd.DataFrame()
    for name,group in bondLenghtsTrajectory(lammpsData, trajectory, bondType).groupby('bType'):
        print('*',group.dist.reset_index(drop=True))
        result[name] = group.dist.reset_index(drop=True)
        print(result)
        #result.append((name, group.dist.hist()))
    result.hist()
    return result

def bondTypeHist2(lammpsData, trajectory, bondType=None):
    '''Return a list of tuples, with the selected bondtype and the axes of 
    the histogram '''
    result = []
    for name,group in bondLenghtsTrajectory(lammpsData, trajectory, bondType).groupby('bType'):
        result.append((name, group.dist.hist()))
        
    return result

if __name__ == "__main__":
    
    #Prepares the data        
    ld = LammpsData(atoms=AtomsFull(), topology=MolecularTopologyData(), forceField=CharmmForceField(), region=Box())
    namdChignolin = NAMDdata()
    namdChignolin.readFiles('../tests/chignolin/chignolin.pdb', '../tests/chignolin/chignolin.psf', '../tests/chignolin/chignolin.prm')
    
    ld.loadNAMDdata(namdChignolin)
    trayectoria = dcdReader.Trajectory("../tests/chignolin/chignolin.dcd")
    
    bondLgth = ld.bondLength()
    
    
    #print("Output de la tabla para el histograma:\n",describeBondTrajectory(ld, trayectoria,5))
    #print("Output del histograma:\n",bondLenghtsTrajectory(ld, trayectoria,[5,6]))
    #trayectoria = dcdReader.Trajectory("../tests/chignolin/chignolin.dcd")
    #print("Output del histograma:\n",bondTypeHist(ld, trayectoria,5))
    #print(bondLenghtsTrajectory(ld, trayectoria,5).groupby("bType"))
    #print("Output del histograma:\n",bondTypeHist(ld, trayectoria,5))
    #fig = plt.subplots((1,5))
    result = bondLenghtsTrajectory(ld, trayectoria,[5,6])
    #r=bondTypeHist(ld, trayectoria,[5])
    
    
    #print("Output de la tabla para el histograma:\n",describeAngleTrajectory(ld, trayectoria))
    
    
    
    

   





