
'''
 chignolin.py

Tests conversion of nanotube simulation from NAMD to LAMMPS

Lyxaira Glass 2019


Modified January 2020 by
  José O. Sotero Esteva
  (jose.sotero@upr.edu)
  Deartment of Mathematics
  University of Puerto Rico at Humacao
'''


from granules.structure.NAMDdata import PDB, PSF, PRM, NAMDdata
from granules.structure.LAMMPSdata import LammpsData

l = LammpsData()
ch = NAMDdata()
ch.readFiles("tubos.pdb", "tubos.psf", "tubos.prm")

l.loadNAMDdata(ch)
l.writeConf("tubos.data")

#print(l.charmmForce())
