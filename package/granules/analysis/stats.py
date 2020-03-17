#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 15:51:33 2019

@author: jse
"""
import numpy as np
import pandas as pd

class RDF:
    def __init__(self, center, maxdist, bins):
        self.center = center
        self.maxdist = maxdist
        self.bins = pd.DataFrame(np.zeros((bins,)), columns=['count'])
        self.bins.index.name = 'levels'
    
    
    def addFrame(self, atoms):
        # generate all pairs of atoms IDs
        atomIds = atoms[['aID']].copy()
        atomIds['key'] = np.ones(len(atomIds))
        atomIds = pd.merge(atomIds, atomIds, on='key')[['aID_x', 'aID_y']]
        atomIds = atomIds[atomIds['aID_x'] < atomIds['aID_y']]

        # compute pairwise distances
        print(len(atomIds))
        from scipy.spatial.distance import pdist
        atomIds['rij'] = pdist(atoms.set_index('aID')[['x', 'y', 'z']].values)
        atomIds = atomIds[atomIds['rij'] < self.maxdist] # remove distant atoms

        #compute distances 
        dists = atomIds['rij']
        #dists = np.sqrt(np.sum(
        #        (atoms[['x','y','z']] - self.center)** 2, 
        #        axis='columns'))
        #dists = dists[dists < self.maxdist]  # remove distant atoms
        delta = self.maxdist / len(self.bins)
        #levels = (dists/(2 * np.pi * delta**3 * len(atoms)**2)).astype('int32').rename('levels').to_frame()
        levels = (dists/delta).astype('int32').rename('levels').to_frame()
        levels = levels.join(pd.DataFrame({'count':np.ones((len(levels),))}))

        r = (levels.levels + 0.5) * delta
        print(r**2)
        #print(levels.groupby('levels').sum() / len(atoms)**2 * 4 * np.pi * r**2)
        #newBins = levels.groupby('levels').sum() * 4 * np.pi * r**2 * delta / len(atoms)**2
        
        
        #self.bins = self.bins.add(newBins, fill_value=0)
       
    def addFrame2(self, atoms):
        #compute distances to center
        dists = np.sqrt(np.sum(
                (atoms[['x','y','z']] - self.center)** 2, 
                axis='columns'))
        dists = dists[dists < self.maxdist]  # remove distant atoms
        delta = self.maxdist / len(self.bins)
        levels = (dists/delta).astype('int32').rename('levels').to_frame()
        levels = levels.join(pd.DataFrame({'count':np.ones((len(levels),))}))

        newBins = levels.groupby('levels').sum()
        self.bins = self.bins.add(newBins, fill_value=0)
       
    def getHist(self):
        print(self.bins)
        delta = self.maxdist / len(self.bins)
        hist = self.bins.join(pd.DataFrame(self.bins.index * delta))
        hist.rename(columns={'levels':'dists'}, inplace=True)
        return hist

if __name__ == "__main__":  # tests
    from granules.structure.NAMDdata import PDB, PSF, PRM, NAMDdata
    from granules.structure.LAMMPSdata import LammpsData

    l = LammpsData()
    ch = NAMDdata()
    ch.readFiles("../tests/tubes/tubos.pdb", "../tests/tubes/tubos.psf", "../tests/tubes/tubos.prm")
    
    l.loadNAMDdata(ch)
    
    rdf = RDF(np.array([0,0,50]), 50, 50)
    rdf.addFrame(l.atoms)
    h = rdf.getHist()
    
    p = h.plot(y='count', x='dists', fontsize=20)
    p.set_title('Radial Distribution Function', fontsize=20)
    p.set_xlabel('distance', fontsize=20)
    p.set_ylabel('count', fontsize=20)
