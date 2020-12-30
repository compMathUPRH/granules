
'''
 chignolin.py

Tests cnversion from NAMD to LAMMPS

Lyxaira Glass 2019
'''


from granules.structure.NAMDdata import PDB, PSF, PRM, NAMDdata
from granules.structure.LAMMPSdata import LammpsData, AtomsFull, MolecularTopologyData, CharmmForceField, Box


l = LammpsData(atoms=AtomsFull(), topology=MolecularTopologyData(), forceField=CharmmForceField(), region=Box())
ch = NAMDdata()
ch.readFiles("2rvd_autopsf.pdb", "2rvd_autopsf.psf", "par_all36_prot.prm")

l.loadNAMDdata(ch)
l.writeConf("2rvd.data")

print("La distancia de los atomos son: ")
#print(l.atomproperty.atoms)
print("###########Bonds###########")
print(l.bondLength())
print(l.bondLength(bondsID=l.topology.bonds.selectTypes(2)))
print("\n",l.meanBondLengthByType())
print(l.describeBonds())
print("###########Angles###########")
print(l.angleLength(l.topology.angles.anID.values))
print(l.meanAngleLengthByType())
print(l.describeAngle())
#print(l.topology.bonds)

'''
print("Tabla de coordenadas: ")
print(l.atomproperty.atoms[['x','y','z']])
print("Tabla enlaces: ")
print(l.topology.bonds)


#Pruebas para selectAtoms()

print("Tabla antes de buscar el atomo: #######")
#print(l.atomproperty.atoms)
print("#####################bonds")
#print(l.topology.bonds)
#print(len(l.topology.bonds))
print("#####################angles")
print(l.topology.angles)
print(len(l.topology.angles))
print("#####################impropers")
#print(l.topology.impropers)
#print(len(l.topology.impropers))
print("#####################dihedrals")
#print(l.topology.dihedrals)
#print(l.forceField.angleCoeffs)
#print(len(l.forceField.angleCoeffs))
l.selectAtom(2)
print("Despues########")
print(l.topology.angles)
print(len(l.topology.angles))
#print(l.topology.dihedrals.head())
#print(len(l.topology.dihedrals))
'''


'''
##Prueba adicional: improper,dihedral
#improper
print(l.topology.impropers)
print(len(l.topology.impropers))

l.selectAtom(2)

#dihedral
print(l.topology.dihedrals)
print(len(l.topology.dihedrals))

##Prueba para atoms y bonds
#print(l.atomproperty.atoms.head())
print("Tabla antes de buscar el atomo: #######")
print(l.topology.bonds)
print(len(l.topology.bonds))
l.selectAtom(2)
#print(l.atomproperty.atoms.tail())
print("Despues########")
print(l.topology.bonds)
print(len(l.topology.bonds))

#Prueba para calcular center() y centerMass()
<<<<<<< HEAD
#print(l.charmmForce())
=======
print("l.charmmForce() = ", l.charmmForce())
c = l.atomproperty.atoms.center()
print(c)
print(l.atomproperty.atoms)
l.atomproperty.atoms.move(-c)
print(l.atomproperty.atoms)
>>>>>>> 6e625574efc27189cd39e6a059963996669764f8
'''
