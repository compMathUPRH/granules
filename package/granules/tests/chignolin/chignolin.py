
'''
 chignolin.py

Tests cnversion from NAMD to LAMMPS

Lyxaira Glass 2019
'''


from granules.structure.NAMDdata import PDB, PSF, PRM, NAMDdata
from granules.structure.LAMMPSdata import LammpsData

l = LammpsData()
ch = NAMDdata()
ch.readFiles("2rvd_autopsf.pdb", "2rvd_autopsf.psf", "par_all36_prot.prm")

l.loadNAMDdata(ch)
l.writeConf("2rvd.data")

print(l.charmmForce())
c = l.atomproperty.atoms.center()
print(c)
print(l.atomproperty.atoms)
l.atomproperty.atoms.move(-c)
print(l.atomproperty.atoms)

