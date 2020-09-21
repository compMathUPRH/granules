
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
import os

l = LammpsData()
ch = NAMDdata()
ch.readFiles("tubos.pdb", "tubos.psf", "tubos.prm")

l.loadNAMDdata(ch)
l.writeConf("tubos.data")

print(l.charmmForce())
#subprocess.call(["lammps", "-in", "tubos.data"])#######Error (arreglar)

'''codigo prueba 
#print(l.charmmForce())
#subprocess.call(["lammps", "-in", "tubos.data"])
print(len(l.atomproperty.atoms))
l.atomproperty.atoms.updateCoordinates("dump.tube")
print(len(l.atomproperty.atoms))

#prueba anterior(vieja)
print(l.atomproperty.atoms[l.atomproperty.atoms['aID'] == 8197])
os.system("lammps -in in.tubos")
l.atomproperty.atoms.updateCoordinates("dump.tube")
print(l.atomproperty.atoms[l.atomproperty.atoms['aID'] == 8197])
'''
