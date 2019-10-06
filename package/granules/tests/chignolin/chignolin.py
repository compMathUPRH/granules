
'''
 chignolin.py

Tests cnversion from NAMD to LAMMPS

Lyxaira Glass 2019
'''


from NAMDdata import NAMDdata
from LAMMPSdata import LammpsData

l = LammpsData()
ch = NAMDdata()
ch.readFiles("2rvd_autopsf.pdb", "2rvd_autopsf.psf", "par_all36_prot.prm")

l.loadNAMDdata(ch)
l.writeConf("2rvd.data")
