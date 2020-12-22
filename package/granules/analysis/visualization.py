#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 21:32:03 2020

@author: lemuel
"""

#heatMap

from granules.simulation import dcdReader
from granules.structure.LAMMPSdata import LammpsData, AtomsFull, MolecularTopologyData, CharmmForceField, Box
from granules.structure.NAMDdata import NAMDdata
import pandas as pd
import matplotlib.pyplot as plt


namdChignolin = NAMDdata()
namdChignolin.readFiles('../tests/chignolin/chignolin.pdb', '../tests/chignolin/chignolin.psf', '../tests/chignolin/chignolin.prm')
#namdChignolin = namdChignolin.pdb[namdChignolin.pdb.ResName != 'WAT']

#select backbone
nombre_atomos = ['CA', 'C', 'N']
namdChignolin.select({'Name':nombre_atomos})

ld = LammpsData(atoms=AtomsFull(), topology=MolecularTopologyData(), forceField=CharmmForceField(), region=Box())

ld.loadNAMDdata(namdChignolin)
trayectoria = dcdReader.Trajectory("../tests/chignolin/chignolin.dcd")

contador = 0
corrida = pd.DataFrame()

print("Número de frame: ", end='')
for X,Y,Z,pbc in trayectoria:
    print("{}, ".format(contador), end='')
    ld.atomproperty.atoms.update_corrdinates((X,Y,Z))
    bonds = ld.bondLength(ld.topologia.bonds.bID.values)  # todos los enlaces
    angles = ld.angleLength(ld.topologia.angles.anID.values)
    corrida[contador] = bonds.dist.append(angles).reset_index(drop=True)
    contador += 1

correlacion = corrida.corr()
temporal = correlacion.replace(1,correlacion.mean())
espacial = corrida.transpose().corr().replace(1,0)

fig, ax = plt.subplots(1,2, figsize=(16,6))
im0 = ax[0].imshow(temporal, cmap ="coolwarm")
im1 = ax[1].imshow(espacial, cmap ="coolwarm")
ax[0].set_title("Temporal")
ax[1].set_title("Espacial")
cbar = fig.colorbar(im1, ax=ax.ravel().tolist(), shrink=0.55)
cbar = fig.colorbar(im0, ax=ax.ravel().tolist(), shrink=0.55)
