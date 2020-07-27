# -*- coding: utf-8 -*-
#---------------------------------------------------------------------------
"""
  LAMMPSdata.py
  First version, October, 2011


  Part of granules Version 0.1.0, October, 2019
    Copyright 2019: José O.  Sotero Esteva, Lyxaira M. Glass Rivera
    Computational Science Group, Department of Mathematics, 
    University of Puerto Rico at Humacao 
    <jose.sotero@upr.edu>.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License version 3 as published by
    the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program (gpl.txt).  If not, see <http://www.gnu.org/licenses/>.

    Acknowledgements: The main funding source for this project has been provided
    by the UPR-Penn Partnership for Research and Education in Materials program, 
    USA National Science Foundation grant number DMR-0934195. 
"""



import pandas as pd
import numpy as np
import random

try:
    import networkx as nx
except Exception as e:
    print("module networkx not found, some functionality may not be available")

def findWithX(typeTuple, pardict):
    ''' Pair 'typeTuple' tuples with correct coefficients found on the 'pardict' 
        dictionary. If tuple is not found, verify dictionary for 'X' values 
        because they can substitute any other type.        
    '''
    try:
        return pardict[typeTuple]
    except KeyError:
        #print("except")
        for k in pardict.keys():
            match = True
            for i in range(len(typeTuple)):
                match = match and (k[i] == 'X' or typeTuple[i] == k[i])
            #print(k,typeTuple, match)
            if match: return pardict[k]
        return np.nan



def detectAtomTypes(charmm):
    ''' extracts atom types from charm.psf.atoms and buils a charmm-lammps
        type translation dictionary.

    Parameter
    -------------
    charmm : NAMDdata
        a NAMDdata object

    Returns:
        a translation dictionary. keys are charmm types, values are integers.
    '''

    # detect different atom types 
    atom_types      = charmm.psf.atoms[['ID','Type']].copy()
    types           = atom_types.drop_duplicates(subset='Type')
    types           = pd.DataFrame(data={'Type':types['Type'].values}) #reindexes
    types['ID']     = types.index + 1
    types           .set_index('Type', inplace=True)
    return types.to_dict()['ID']


class LammpsBodySection(pd.DataFrame):

    def add(self, data):
    
        #Almacenar numeros de columnas
        columns = list(self.columns)  
        
        #Separar la data para almacenar en diccionario
        dtable = {col:[] for col in columns}
        
        #Rellenar lista de listas con valor nulo '0' de ser necesario
        dif = len(self.columns) - len(data[0])
        if dif > 0: 
            for i in range(len(data)): 
                for j in range(dif):
                    data[i].append('0')
        
        for d in data:
            for col in columns:
                try:dtable[col].append(self[col].dtype.type(d.pop(0)))
                except KeyError:dtable[col].append(float(d.pop(0)))

        #Crear dataframe con toda la informacion
        fr = pd.DataFrame(dtable, copy=False, columns=columns)
       
        #Ajustar indices para que no comiencen en 0
        fr.index = np.arange(1, len(fr) + 1)
       
        #Reescribir DataFrame original con creado
        if self.empty:
            super().__init__(data=fr)
        else:
            super().__init__(data=self.append(fr, ignore_index=True))


#===================================================================

class AtomPropertyData():
    def __init__(self):
        
        # atom-property sections
        self.atoms       = AtomsDF()
        self.velocities  = VelocitiesDF()
        self.masses      = MassesDF()
               
    def setFromNAMD(self,charmm): 
        '''Llama a la funcion setFromNAMD() de las clases de la clase AtomPropertyData,
            asignadas en los atributo.'''
            
        self.atoms.setFromNAMD(charmm)
        self.velocities.setToZero(self.atoms)
        self.masses.setFromNAMD(charmm.psf.atoms,self.atoms)

    def centerMass(self):
        '''Regresa el centro de masa.'''
        coordMass = self.atoms.set_index('aType')[['x','y','z']].join(self.masses.set_index('aType'))
        for c in ['x','y','z']:
            coordMass[c] = coordMass[c] * coordMass['Mass']
        coordMass = coordMass[['x','y','z']]

        return coordMass.mean()
    
   
  
class MolecularTopologyData():
    def __init__(self):
        
        # molecular topology sections
        self.angles      = AnglesDF()
        self.bonds       = BondsDF()
        self.dihedrals   = DihedralsDF()
        self.impropers   = ImpropersDF()
        
    def setFromNAMD(self,charmm): 
        '''Llama a la funcion setFromNAMD() de las clases de la clase MolecularTopolyData,
            asignadas en los atributo.'''
        
        # molecular topology sections
        self.bonds.setFromNAMD(charmm)
        self.angles.setFromNAMD(charmm)
        self.dihedrals.setFromNAMD(charmm)
        self.impropers.setFromNAMD(charmm)

  

class ForceFieldData():
     def __init__(self):

        #force field sections
        #self.angleAngleCoeffs       
        #self.angleAngleTorsionCoeffs
        #self.angleTorsionCoeffs     
	
        #self.bondAngleCoeffs        
        #self.bondBond13Coeffs       
        
        #self.middleBondTorsion      
        #self.endBondTorsionCoeffs  
        self.angleCoeffs           = AngleCoeffs()
        self.bondCoeffs             = BondCoeffs()
        
        self.dihedralCoeffs         = DihedralCoeffs()
        self.improperCoeffs         = ImproperCoeffs()
        self.pairCoeffs             = PairCoeffs()
        
     def setFromNAMD(self,charmm,atompropertydata,topology): #añadi este codigo nuevo 
        
        self.pairCoeffs.setFromNAMD(charmm, atompropertydata.masses)
        self.bondCoeffs.setFromNAMD(charmm, topology.bonds)
        self.angleCoeffs.setFromNAMD(charmm, topology.angles)
        self.dihedralCoeffs.setFromNAMD(charmm, topology.dihedrals)
        self.improperCoeffs.setFromNAMD(charmm, topology.impropers)
        
     def charmmNonBondEnergy(self,atompropertydata,topologia):
        ''' Computes CHARMM Lennard-Jones energy.
            Formula: Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
                    Eps,i,j = sqrt(eps,i * eps,j)
                    Rmin,i,j = Rmin/2,i + Rmin/2,j

            Computes CHARMM Coulumb energy.
            Formula: qi qj/rij / epsilon_0

            returns (L-J, Coulomb)
        '''
        from scipy.constants import epsilon_0, physical_constants

        NONB_CUTOFF = 13.0

        print("ForceFieldData.charmmNonBondEnergy")
        
        # generate all pairs of atoms IDs
        atoms = atompropertydata.atoms.copy() #Cambio
        atomIds = atoms[['aID']]
        print("ForceFieldData.charmmNonBondEnergy atomIds.columns=", atomIds.columns)
        atomIds = atomIds.assign(key=np.ones(len(atomIds)))
        #atomIds['key'] = np.ones(len(atomIds))
        #print("ForceFieldData.charmmNonBondEnergy atomIds.columns=", atomIds.columns)
        atomIds = pd.merge(atomIds, atomIds, on='key')[['aID_x', 'aID_y']]
        atomIds = atomIds[atomIds['aID_x'] < atomIds['aID_y']]

        atomIds['nbID'] = np.arange(len(atomIds))

        # compute pairwise distances
        print("ForceFieldData.charmmNonBondEnergy len(atomIds) =", len(atomIds))
        from scipy.spatial.distance import pdist
        atomIds['rij'] = pdist(atoms.set_index('aID')[['x', 'y', 'z']].values)
        atomIds = atomIds[atomIds['rij'] < NONB_CUTOFF]

        # remove bonded atoms
        atomIds['p'] = list(zip(atomIds.aID_x, atomIds.aID_y))
        bonds = topologia.bonds.copy()
        bonds['p'] = list(zip(bonds.Atom1, bonds.Atom2))
        atomIds = atomIds.set_index('p').join(bonds.set_index('p'))
        atomIds = atomIds[atomIds.bID.isna()][['aID_x','aID_y','nbID','rij']].reset_index(drop=True)
        del bonds

        # remove angled atoms
        atomIds['p'] = list(zip(atomIds.aID_x, atomIds.aID_y))
        angles = topologia.angles.copy()
        angles['p'] = list(zip(angles.Atom1, angles.Atom3))
        atomIds = atomIds.set_index('p').join(angles.set_index('p'))
        atomIds = atomIds[atomIds.anID.isna()][['aID_x','aID_y','nbID','rij']].reset_index(drop=True)
        del angles

        # get atom types and charges
        atomIds = atomIds.set_index('aID_x').join(atompropertydata.atoms[['aID', 'Q']].set_index('aID'))
        atomIds.rename(columns={'Q':'qi'}, inplace=True)
        atomIds = atomIds.join(atompropertydata.atoms[['aID', 'aType']].set_index('aID')).reset_index(drop=True)
        atomIds.rename(columns={'aType':'aiType'}, inplace=True)
        atomIds = atomIds.set_index('aID_y').join(atompropertydata.atoms[['aID', 'Q']].set_index('aID'))
        atomIds = atomIds.join(atompropertydata.atoms[['aID', 'aType']].set_index('aID')).reset_index(drop=True)
        atomIds.rename(columns={'aType':'ajType', 'Q':'qj'}, inplace=True)
        print("ForceFieldData.charmmNonBondEnergy len(atomIds (clean)) =", len(atomIds))

        # get epsilons and sigmas for each atom type
        atomIds = atomIds.set_index('aiType').join(
                        self.pairCoeffs.set_index('aType')
                 ).reset_index(drop=True)
        atomIds.drop(columns=['epsilon1_4', 'sigma1_4'], inplace=True)
        atomIds.rename(columns={'epsilon':'epsilon_i', 'sigma':'sigma_i'}, inplace=True)
        atomIds = atomIds.set_index('ajType').join(
                        self.pairCoeffs.set_index('aType')
                 ).reset_index(drop=True).drop(columns=['epsilon1_4', 'sigma1_4'])
        atomIds.rename(columns={'epsilon':'epsilon_j', 'sigma':'sigma_j'}, inplace=True)

        # compute epsilon and sigma
        atomIds['epsilon'] = np.sqrt(atomIds.epsilon_i * atomIds.epsilon_j)
        atomIds['sigma'] = atomIds.sigma_i + atomIds.sigma_j
        atomIds.drop(columns=['epsilon_i', 'epsilon_j', 'sigma_i', 'sigma_j'], inplace=True)


        atomIds.set_index('nbID', inplace=True)
        print("ForceFieldData.charmmNonBondEnergy END")
        #print(atomIds)
        # return LJ and Coulomb
        COULOMB = 332.0636

        '''
        return -np.sum(atomIds.epsilon * (atomIds.rij**12 - atomIds.rij**6)), \
               np.sum((physical_constants['electric constant'][0] * atomIds.qi * atomIds.qj) / (atomIds.rij * epsilon_0))
        '''
        return -np.sum(atomIds.epsilon * ((atomIds.sigma/atomIds.rij)**12 - (atomIds.sigma/atomIds.rij)**6)), \
               COULOMB * np.sum((atomIds.qi * atomIds.qj) / (atomIds.rij))

     def charmmBondEnergy(self,atompropertydata,topologia):
        ''' Computes CHARMM bond energy.

            Formula: sum K * (bij - b0)**2
        '''
        bi = topologia.bonds.set_index('Atom1').join(
                atompropertydata.atoms.set_index('aID')
             )[['bID','x', 'y', 'z']].set_index('bID')

        bj = topologia.bonds.set_index('Atom2').join(
                atompropertydata.atoms.set_index('aID')
             )[['bID','x', 'y', 'z']].set_index('bID')


        bij = np.sqrt(np.sum((bi - bj) ** 2, axis='columns'))  # bonds lengths

        coeffs = topologia.bonds[['bID','bType']].set_index('bType').join(
                        self.bondCoeffs.set_index('bType')
                 ).reset_index(drop=True)
        K = coeffs[['bID','Spring_Constant']].set_index('bID').Spring_Constant
        b0 = coeffs[['bID','Eq_Length']].set_index('bID').Eq_Length

        return np.sum(K * (bij-b0)**2)


     def charmmAngleEnergy(self,atompropertydata,topologia):
        ''' Computes CHARMM angle energy.
            Formula: sum K * (aij - a0)**2

        '''
        bi = topologia.angles.set_index('Atom1').join(
                atompropertydata.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        bj = topologia.angles.set_index('Atom2').join(
                atompropertydata.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        bk = topologia.angles.set_index('Atom3').join(
                atompropertydata.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        # compute angles
        l1 = bi - bj
        l2 = bk - bj

        norm1 = np.sqrt(np.square(l1).sum(axis=1))
        norm2 = np.sqrt(np.square(l1).sum(axis=1))
        dot = l1.x * l2.x + l1.y * l2.y + l1.z * l2.z
        angles = np.arccos(dot / (norm1 * norm2))

        coeffs = topologia.angles[['anID','anType']].set_index('anType').join(
                        self.angleCoeffs.set_index('anType')
                 ).reset_index(drop=True)

        K = coeffs[['anID','Ktheta']].set_index('anID').Ktheta
        a0 = np.radians(coeffs[['anID','Theta0']].set_index('anID').Theta0)

        return np.sum(K * (angles-a0)**2)

     def charmmDihedralsEnergy(self,atompropertydata,topologia):
        ''' Computes CHARMM angle energy.
            Formula: sum K * (1 + cos(n * x - d))
        '''
        bi = topologia.dihedrals.set_index('Atom1').join(
                atompropertydata.atoms.set_index('aID')
             )[['dID','x', 'y', 'z']].set_index('dID')

        bj = topologia.dihedrals.set_index('Atom2').join(
                atompropertydata.atoms.set_index('aID')
             )[['dID','x', 'y', 'z']].set_index('dID')

        bk = topologia.dihedrals.set_index('Atom3').join(
                atompropertydata.atoms.set_index('aID')
             )[['dID','x', 'y', 'z']].set_index('dID')

        bl = topologia.dihedrals.set_index('Atom4').join(
                atompropertydata.atoms.set_index('aID')
             )[['dID','x', 'y', 'z']].set_index('dID')

        # compute equations of planes
        l1 = bi - bj
        l2 = bk - bj
        normals1 = np.cross(l1, l2)
        normals1 = normals1 / np.linalg.norm(normals1, axis=1)[:, np.newaxis] 

        l1 = bj - bk
        l2 = bl - bk
        normals2 = np.cross(l1, l2)
        normals2 = normals2 / np.linalg.norm(normals2, axis=1)[:, np.newaxis] 

        angles = np.degrees(np.arccos(np.abs(np.einsum('ij,ij->i', normals1, normals2))))

        # get CHARMM constants
        print(topologia.dihedrals)
        coeffs = topologia.dihedrals[['dID','dType']].set_index('dType').join(
                        self.dihedralCoeffs.set_index('dType')
                 ).reset_index(drop=True)
        K = coeffs[['dID','Kchi']].set_index('dID').Kchi
        n = coeffs[['dID','n']].set_index('dID').n
        d = coeffs[['dID','delta']].set_index('dID').delta



        return np.sum(K * (1 + np.cos(n * angles - d)))
        '''
        '''


     def charmmNonBondForce(self,atompropertydata,topologia):
        ''' Computes CHARMM Lennard-Jones energy.
            Formula: Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
                    Eps,i,j = sqrt(eps,i * eps,j)
                    Rmin,i,j = Rmin/2,i + Rmin/2,j

            Computes CHARMM Coulumb energy.
            Formula: qi qj/rij / epsilon_0

            returns (L-J, Coulomb)
        '''
        from scipy.constants import epsilon_0, physical_constants

        NONB_CUTOFF = 13.0
        print("ForceFieldData.charmmNonBondForce()")
        
        # generate all pairs of atoms IDs
        atomIds = atompropertydata.atoms[['aID']].copy()
        atomIds['key'] = np.ones(len(atomIds))
        atomIds = pd.merge(atomIds, atomIds, on='key')[['aID_x', 'aID_y']]
        atomIds = atomIds[atomIds['aID_x'] < atomIds['aID_y']]

        atomIds['nbID'] = np.arange(len(atomIds))

        # compute pairwise distances
        print('ForceFieldData.charmmNonBondForce: len(atomIds)=', len(atomIds))
        from scipy.spatial.distance import pdist
        atomIds['rij'] = pdist(atompropertydata.atoms.set_index('aID')[['x', 'y', 'z']].values)
        atomIds = atomIds[atomIds['rij'] < NONB_CUTOFF]

        # remove bonded atoms
        print('ForceFieldData.charmmNonBondForce: len(atomIds < NONB_CUTOFF)=', len(atomIds))
        atomIds['p'] = list(zip(atomIds.aID_x, atomIds.aID_y))
        bonds = topologia.bonds.copy()
        bonds['p'] = list(zip(bonds.Atom1, bonds.Atom2))
        atomIds = atomIds.set_index('p').join(bonds.set_index('p'))
        atomIds = atomIds[atomIds.bID.isna()][['aID_x','aID_y','nbID','rij']].reset_index(drop=True)
        del bonds
        print('ForceFieldData.charmmNonBondForce: len(atomIds < NONB_CUTOFF -  BONDS)=', len(atomIds))

        # remove angled atoms
        atomIds['p'] = list(zip(atomIds.aID_x, atomIds.aID_y))
        angles = topologia.angles.copy()
        angles['p'] = list(zip(angles.Atom1, angles.Atom3))
        atomIds = atomIds.set_index('p').join(angles.set_index('p'))
        atomIds = atomIds[atomIds.anID.isna()][['aID_x','aID_y','nbID','rij']].reset_index(drop=True)
        del angles
        print('ForceFieldData.charmmNonBondForce: len(atomIds < NONB_CUTOFF -  BONDS - ABGLES)=',len(atomIds))

        # get atom types and charges
        atomIds = atomIds.set_index('aID_x', drop=False).join(atompropertydata.atoms[['aID', 'Q']].set_index('aID'))
        atomIds.rename(columns={'Q':'qi'}, inplace=True)
        atomIds = atomIds.join(atompropertydata.atoms[['aID', 'aType']].set_index('aID')).reset_index(drop=True)
        atomIds.rename(columns={'aType':'aiType'}, inplace=True)
        atomIds = atomIds.set_index('aID_y', drop=False).join(atompropertydata.atoms[['aID', 'Q']].set_index('aID'))
        atomIds = atomIds.join(atompropertydata.atoms[['aID', 'aType']].set_index('aID')).reset_index(drop=True)
        atomIds.rename(columns={'aType':'ajType', 'Q':'qj'}, inplace=True)

        # get epsilons and sigmas for each atom type
        atomIds = atomIds.set_index('aiType').join(
                        self.pairCoeffs.set_index('aType')
                 ).reset_index(drop=True)
        atomIds.drop(columns=['epsilon1_4', 'sigma1_4'], inplace=True)
        atomIds.rename(columns={'epsilon':'epsilon_i', 'sigma':'sigma_i'}, inplace=True)
        atomIds = atomIds.set_index('ajType').join(
                        self.pairCoeffs.set_index('aType')
                 ).reset_index(drop=True).drop(columns=['epsilon1_4', 'sigma1_4'])
        atomIds.rename(columns={'epsilon':'epsilon_j', 'sigma':'sigma_j'}, inplace=True)

        # compute epsilon and sigma
        atomIds['epsilon'] = np.sqrt(atomIds.epsilon_i * atomIds.epsilon_j)
        atomIds['sigma'] = atomIds.sigma_i + atomIds.sigma_j
        atomIds.drop(columns=['epsilon_i', 'epsilon_j', 'sigma_i', 'sigma_j'], inplace=True)


        atomIds.set_index('nbID', drop=False, inplace=True)
        # return LJ and Coulomb
        COULOMB = 332.0636


        bi = atomIds.set_index('aID_x').join(
                atompropertydata.atoms.set_index('aID')
             )[['nbID', 'x', 'y', 'z']].set_index('nbID')
        bj = atomIds.set_index('aID_y').join(
                atompropertydata.atoms.set_index('aID')
             )[['nbID', 'x', 'y', 'z']].set_index('nbID')
        
        bij = (bj - bi).div(atomIds.rij, axis=0)
        
        forces = -atomIds.epsilon * (-12 * (atomIds.sigma**12/atomIds.rij**13) + 6 * (atomIds.sigma**6/atomIds.rij**7)) - COULOMB * (atomIds.qi * atomIds.qj) / (atomIds.rij**2)
        wi = bij.mul(forces, axis=0)
        wj = bij.mul(-forces, axis=0)
        
        fi = wi.join(atomIds[['nbID', 'aID_x']].set_index('nbID')).groupby('aID_x').sum()
        fj = wj.join(atomIds[['nbID', 'aID_y']].set_index('nbID')).groupby('aID_y').sum()
        ff = fj.add(fi, axis=0, fill_value=0)

        return ff


     def charmmBondForce(self,atompropertydata,topologia):
        ''' Computes CHARMM bond energy.
            Formula: sum K * (bij - b0)**2

        '''
        bi = topologia.bonds.set_index('Atom1').join(
                atompropertydata.atoms.set_index('aID')
             )[['bID','x', 'y', 'z']].set_index('bID')

        bj = topologia.bonds.set_index('Atom2').join(
                atompropertydata.atoms.set_index('aID')
             )[['bID','x', 'y', 'z']].set_index('bID')

        bij = bj - bi  # bonds
        rbij = np.sqrt(np.sum((bij) ** 2, axis='columns'))  # bonds lengths

        coeffs = topologia.bonds[['bID','bType']].set_index('bType').join(
                        self.bondCoeffs.set_index('bType')
                 ).reset_index(drop=True)
        K = coeffs[['bID','Spring_Constant']].set_index('bID').Spring_Constant
        b0 = coeffs[['bID','Eq_Length']].set_index('bID').Eq_Length

        forces = -2 * K * (rbij-b0)
        
        bij = bij.div(rbij, axis=0)  # normalize
        wi = bij.mul(forces, axis=0)
        wj = bij.mul(-forces, axis=0)
        
        fi = wi.join(topologia.bonds.set_index('bID'))[['x','y','z','Atom1']].groupby('Atom1').sum()
        fj = wi.join(topologia.bonds.set_index('bID'))[['x','y','z','Atom2']].groupby('Atom2').sum()
        fi.index.names = ['aID']
        fj.index.names = ['aID']
        ff = fi.add(fj, axis=0, fill_value=0)

        return ff



     def charmmAngleForce(self,atompropertydata,topologia):
        ''' Computes CHARMM angle energy.
            Formula: sum K * (aij - a0)**2
        '''
        
        print("ForceFieldData.charmmAngleForce")
        bi = topologia.angles.set_index('Atom1').join(
                atompropertydata.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        bj = topologia.angles.set_index('Atom2').join(
                atompropertydata.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        bk = topologia.angles.set_index('Atom3').join(
                atompropertydata.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        # compute angles
        l1 = bi - bj
        l2 = bk - bj

        norm1 = np.sqrt(np.square(l1).sum(axis=1))
        norm2 = np.sqrt(np.square(l1).sum(axis=1))
        dot = l1.x * l2.x + l1.y * l2.y + l1.z * l2.z
        angles = np.arccos(dot / (norm1 * norm2))

        coeffs = topologia.angles[['anID','anType']].set_index('anType').join(
                        self.angleCoeffs.set_index('anType')
                 ).reset_index(drop=True)

        K = coeffs[['anID','Ktheta']].set_index('anID').Ktheta
        a0 = np.radians(coeffs[['anID','Theta0']].set_index('anID').Theta0)

        c11 = l1.x * l1.x + l1.y * l1.y + l1.z * l1.z
        c12 = l1.x * l2.x + l1.y * l2.y + l1.z * l2.z
        c22 = l2.x * l2.x + l2.y * l2.y + l2.z * l2.z
        cd = np.sqrt(c11 * c22)
        c = c12 / cd
        
        forces = -2 * K * (angles-a0)
        
        w1 = l1.mul(c12 / c11, axis=0).add(-l2, axis=0).mul(forces / cd, axis=0)
        w2 = l1.sub(l2.mul(c12 / c22, axis=0),  axis=0).mul(forces / cd, axis=0)

        #print(self.angles.columns.values)
        f1 = w1.join(topologia.angles.set_index('anID'))[['x','y','z','Atom1']].groupby('Atom1').sum()
        f2 = w1.sub(w2, axis=0).join(topologia.angles.set_index('anID'))[['x','y','z','Atom2']].groupby('Atom2').sum()
        f3 = w2.join(topologia.angles.set_index('anID'))[['x','y','z','Atom3']].groupby('Atom3').sum()
        f1.index.names = ['aID']
        f2.index.names = ['aID']
        f3.index.names = ['aID']

        ff = f1.add(f2.add(f3, axis=0, fill_value=0), axis=0, fill_value=0)
        #print (ff)
        print("ForceFieldData.charmmAngleForce  END")
        return ff

     def charmmForce(self,atompropertydata,topologia):
        print("ForceFieldData.charmmForce()")
        return self.charmmNonBondForce(atompropertydata,topologia).add(
                self.charmmBondForce(atompropertydata,topologia), axis=0).add(
                self.charmmAngleForce(atompropertydata,topologia), axis=0)
        
     def charmmEnergy(self,atompropertydata,topologia):
        return self.charmmNonBondEnergy(atompropertydata,topologia).add(
                self.charmmBondEnergy(atompropertydata,topologia)).add(
                self.charmmAngleEnergy(atompropertydata,topologia)).add(
                self.charmmDihedralsEnergy(atompropertydata,topologia))

	

#===================================================================
#Clases temporeras que hereda de panda,las subclases de sus clases originales(nombre)
#   se pasaron a estas nuevas clases.

class AtomProperty(LammpsBodySection):#cambiado LammpsBodySection
    pass
'''
    def __init__(self):
        pass     
 '''
   
class MolecularTopology(LammpsBodySection):
    pass
'''
    def __init__(self):
        pass
'''

class ForceField(LammpsBodySection):
    pass
'''
    def __init__(self):
        pass
'''
#=====================================================================
class AtomsDF(AtomProperty):
    def __init__(self,data=None, dtype=None, copy=False):
        dtypes = {'aID':[0], 'Mol_ID':[0], 'aType':[0], 'Q':[0.0], 
                  'x':[0.0], 'y':[0.0], 'z':[0.0], 'Nx':[0], 'Ny':[0], 'Nz':[0]}   
        super(AtomsDF, self).__init__(data=dtypes, copy=copy, columns=['aID', 'Mol_ID', 'aType', 'Q', 
                  'x', 'y', 'z', 'Nx', 'Ny', 'Nz'])
        super(AtomsDF, self).__init__(self.drop([0]))
        
    def setFromNAMD(self, charmm):
        ''' Extracts info from ATOMS object of a PSF object into self.

        Parameter
        -----------------
        charmm : NAMDdata
            NAMDdata object
        '''

        charmmTypeToInt = detectAtomTypes(charmm)
        #print(charmmTypeToInt)
        #print(atom_types['Type'].map(charmmTypeToInt))

        # extract info from charmm
        sel_psf     = charmm.psf.atoms[['ID', 'Type', 'Charge']].set_index('ID')
        #print(sel_psf)
        sel_pdb     = charmm.pdb[['ID','x','y','z']].set_index('ID')
        #print(sel_pdb)
        sel         = sel_pdb.join(sel_psf)
        #print(sel)        
        sel['aType'] = sel_psf['Type'].map(charmmTypeToInt)
        #print(sel_pdb[charmm.pdb['ID']==0])
        sel.reset_index(inplace=True)
        #print(sel)        
        sel         .rename(columns={"Charge":"Q",'ID':'aID'}, inplace=True)

        # add remining columns
        sel['Mol_ID'] = np.ones((len(sel), 1), dtype=np.int8)
        sel['Nx']     = np.zeros((len(sel), 1))
        sel['Ny']     = np.zeros((len(sel), 1))
        sel['Nz']     = np.zeros((len(sel), 1))
        sel['ID']     = np.arange(1, len(sel)+1)
        #print('setFromNAMD:\n', charmmTypeToInt,sel_psf['Type'],sel['aType'])
        #print('setFromNAMD:\n', charmmTypeToInt,sel_psf.loc[1,'Type'],sel_psf.loc[1,'Type'] in charmmTypeToInt, sel.loc[1,'aType'])
        #print('\nsetFromNAMD:\n', sel[np.isnan(sel['aType'])])
        # rearrange columns
        sel = sel[['aID', 'Mol_ID', 'aType', 'Q', 'x', 'y', 'z', 'Nx', 'Ny', 'Nz']]
        
        #sel.reset_index(inplace=True)
        #print("sel = ", sel.dtypes)
        super(AtomsDF, self).__init__(sel.astype({
                     'Mol_ID' :int,
                     'aType' :int,
                     'Q' :float,
                     'x' :float,
                     'y' :float,
                     'z' :float,
                     'Nx' :int,
                     'Ny' :int,
                     'Nz' : int
                    }))
        #print(self.dtypes)
    
    def center(self):
        return self[['x', 'y', 'z']].mean()

    def move(self, b):
        for coord in ['x', 'y', 'z']:
            self[coord] = self[coord] + b[coord]
        return self

    def updateCoordinates(self,archivo):
        '''Actualiza las coordenadas x,y,z del dataframe'''
        fill = open(archivo,"r")
        coor = []
        
        for i in reversed(fill.readlines()):#crea una lista del archivo en reversa
            if 'id' in i: break
            coor.append(i.rstrip("\n").split(" "))
            
        data = pd.DataFrame(coor,columns = ['aID', 'aType', 'x', 'y', 'z'])#crea un dataframe
        data = data.astype({'aID':int, 'aType':int, 'x':float, 'y':float, 'z':float})#type de dato de cada columna
        data = data.sort_values('aID').reset_index(drop=True)
        fill.close()
  
        #reemplaza los valores del dataframe viejo a los valores del dataframe nuevo
        for column in ['x', 'y', 'z']: self[column] = data[column].values
     
class MassesDF(AtomProperty):
    def __init__(self,data=None, dtype=None, copy=False):
        if data  is None:
            dtypes = {'aType':[0], 'Mass':[0.0]}
            super(MassesDF, self).__init__(data=dtypes, copy=copy, columns=dtypes.keys())
            super(MassesDF, self).__init__(self.drop([0]))
 
    def setFromNAMD(self, psf_atoms, atoms):
        ''' Extracts info from ATOMS object of a PSF object into self.

        Parameter
        -----------------
        psf_atoms : PSF.ATOMS
            a PSF.ATOMS object

        atoms     : AtomsDF
            an AtomsDF object
        '''

        # extract info from charmm and LAMMPS
        sel_psf  = psf_atoms[['ID', 'Mass']].set_index('ID')
        sel_self = atoms[['aID', 'aType']].set_index('aID').copy()
        sel      = sel_self.join(sel_psf).drop_duplicates().reset_index()

        # rename columns
        sel      .drop(columns='aID', inplace=True)
        #sel      .rename(columns={"Type":"mID"}, inplace=True)
        #print(sel.dtypes)

        super(MassesDF, self).__init__(sel)


class VelocitiesDF(AtomProperty):
    def __init__(self,data=None, dtype=None, copy=False):
        dtypes = {'vID':[0], 'Vx':[0.0], 'Vy':[0.0], 'Vz':[0.0]}
        super(VelocitiesDF, self).__init__(data=dtypes, copy=copy, columns=dtypes.keys())
        super(VelocitiesDF, self).__init__(self.drop([0]))
   
    def setToZero(self, atoms):
        ''' Sets velocitis of atoms to zero.

        Parameter
        -----------------
        atoms     : AtomsDF
            an AtomsDF object
        '''

        # extract info from LAMMPS
        sel = atoms[['aID']].copy().rename(columns={'aID':'vID'})
        #sel.rename(columns={'aID':'vID'}, inplace=True)
        sel['Vx']     = np.zeros((len(sel), 1))
        sel['Vy']     = np.zeros((len(sel), 1))
        sel['Vz']     = np.zeros((len(sel), 1))
        #print("VelocitiesDF sel = ", sel.dtypes)

        super(VelocitiesDF, self).__init__(sel)

#===================================================================

class AnglesDF(MolecularTopology):
    def __init__(self,data=None, dtype=None, copy=False):
        dtypes = {'anID':[0], 'anType':[0], 'Atom1':[0], 'Atom2':[0], 'Atom3':[0]}
        super(AnglesDF, self).__init__(data=dtypes, copy=copy, columns=dtypes.keys())
        super(AnglesDF, self).__init__(self.drop([0]))

    def setFromNAMD(self, charmm):
        ''' Extracts info from ATOMS object of a PSF object into self.

        Parameter
        -----------------
        charmm : NAMDdata
            NAMDdata object
        '''

        # extract info from charmm
        psf_types   = charmm.psf.atoms[['ID', 'Type']].set_index('ID').to_dict()['Type']
        #print(psf_types)

        # substitute atoms numbers with charmm atom types
        angles     = charmm.psf.angles.copy()
        angles['atom1'] = angles['atom1'].map(psf_types)
        angles['atom2'] = angles['atom2'].map(psf_types)
        angles['atom3'] = angles['atom3'].map(psf_types)
        angles['atuple'] = list(zip(angles.atom1, angles.atom2, angles.atom3))
        angles.drop(columns=['atom1', 'atom2', 'atom3'], inplace=True)
        #print(angles)

        # build translation dict
        btypes = angles.copy().drop_duplicates(inplace=False)
        #print(btypes)
        btypes['ID'] = np.arange(1, len(btypes)+1)
        btypes.set_index('atuple', inplace=True)
        btypeToInt = btypes.to_dict()['ID']
        #print(btypeToInt)
        btypes.reset_index(inplace=True)

        # final table
        angles['ID'] = np.arange(1, len(angles)+1)
        angles['Type'] = angles['atuple'].map(btypeToInt)
        angles.drop(columns=['atuple'], inplace=True)
        angles['Atom1'] = charmm.psf.angles.copy()['atom1']
        angles['Atom2'] = charmm.psf.angles.copy()['atom2']
        angles['Atom3'] = charmm.psf.angles.copy()['atom3']

        angles.rename(columns={'ID':'anID', 'Type':'anType'}, inplace=True)
        #print(angles)

        super(AnglesDF, self).__init__(angles)
        #print(self.dtypes)

class BondsDF(MolecularTopology):
    def __init__(self,data=None, dtype=None, copy=False):
        dtypes = {'bID':[0], 'bType':[0], 'Atom1':[0], 'Atom2':[0]}     
        super(BondsDF, self).__init__(data=dtypes, copy=copy, columns=dtypes.keys())
        super(BondsDF, self).__init__(self.drop([0]))

    def setFromNAMD(self, charmm):
        ''' Extracts info from ATOMS object of a PSF object into self.

        Parameter
        -----------------
        charmm : NAMDdata
            NAMDdata object
        '''

        # extract info from charmm
        psf_types   = charmm.psf.atoms[['ID', 'Type']].set_index('ID').to_dict()['Type']

        # substitute atoms numbers with charmm atom types
        bonds     = charmm.psf.bonds.copy()
        bonds['atom1'] = bonds['atom1'].map(psf_types)
        bonds['atom2'] = bonds['atom2'].map(psf_types)
        bonds['atuple'] = list(zip(bonds.atom1, bonds.atom2))
        bonds.drop(columns=['atom1', 'atom2'], inplace=True)

        # build translation dict
        btypes = bonds.copy().drop_duplicates(inplace=False)
        #print(btypes)
        btypes['ID'] = np.arange(1, len(btypes)+1)
        btypes.set_index('atuple', inplace=True)
        btypeToInt = btypes.to_dict()['ID']
        btypes.reset_index(inplace=True)

        # final table
        bonds['ID'] = np.arange(1, len(bonds)+1)
        bonds['Type'] = bonds['atuple'].map(btypeToInt)
        bonds.drop(columns=['atuple'], inplace=True)
        bonds['Atom1'] = charmm.psf.bonds.copy()['atom1']
        bonds['Atom2'] = charmm.psf.bonds.copy()['atom2']

        bonds.rename(columns={'ID':'bID', 'Type':'bType'}, inplace=True)

        super(BondsDF, self).__init__(bonds)
        

       
class DihedralsDF(MolecularTopology):
    def __init__(self,data=None, dtype=None, copy=False):
        dtypes = {'dID':[0], 'dType':[0], 'atom1':[0], 'atom2':[0], 'atom3':[0], 'atom4':[0]}
        super(DihedralsDF, self).__init__(data=dtypes, copy=copy, columns=dtypes.keys())
        super(DihedralsDF, self).__init__(self.drop([0]))

    def setFromNAMD(self, charmm):
        ''' Extracts info from ATOMS object of a PSF object into self.

        Parameter
        -----------------
        charmm : NAMDdata
            NAMDdata object
        '''

        # extract info from charmm
        psf_types   = charmm.psf.atoms[['ID', 'Type']].set_index('ID').to_dict()['Type']
        #print(psf_types)

        # substitute atoms numbers with charmm atom types
        dihes     = charmm.psf.dihedrals.copy()
        #print(charmm.psf.dihedrals)
        dihes['atom1'] = dihes['atom1'].map(psf_types)
        dihes['atom2'] = dihes['atom2'].map(psf_types)
        dihes['atom3'] = dihes['atom3'].map(psf_types)
        dihes['atom4'] = dihes['atom4'].map(psf_types)
        dihes['atuple'] = list(zip(dihes.atom1, dihes.atom2, dihes.atom3, dihes.atom4))
        dihes.drop(columns=['atom1', 'atom2', 'atom3', 'atom4'], inplace=True)
        #print(dihes)

        # build translation dict
        btypes = dihes.copy().drop_duplicates(inplace=False)
        #print(btypes)
        btypes['ID'] = np.arange(1, len(btypes)+1)
        btypes.set_index('atuple', inplace=True)
        btypeToInt = btypes.to_dict()['ID']
        #print(btypeToInt)
        btypes.reset_index(inplace=True)

        # final table
        dihes['ID'] = np.arange(1, len(dihes)+1)
        dihes['Type'] = dihes['atuple'].map(btypeToInt)
        dihes.drop(columns=['atuple'], inplace=True)
        dihes['Atom1'] = charmm.psf.dihedrals.copy()['atom1']
        dihes['Atom2'] = charmm.psf.dihedrals.copy()['atom2']
        dihes['Atom3'] = charmm.psf.dihedrals.copy()['atom3']
        dihes['Atom4'] = charmm.psf.dihedrals.copy()['atom4']

        dihes.rename(columns={'ID':'dID', 'Type':'dType'}, inplace=True)
        #print(dihes)

        super(DihedralsDF, self).__init__(dihes)
        #print(self.dtypes)


class ImpropersDF(MolecularTopology):
    def __init__(self,data=None, dtype=None, copy=False):
        dtypes = {'iID':[0], 'iType':[0], 'atom1':[0], 'atom2':[0], 'atom3':[0], 'atom4':[0]}
        super(ImpropersDF, self).__init__(data=dtypes, copy=copy, columns=dtypes.keys())
        super(ImpropersDF, self).__init__(self.drop([0]))

    def setFromNAMD(self, charmm):
        ''' Extracts info from ATOMS object of a PSF object into self.

        Parameter
        -----------------
        charmm : NAMDdata
            NAMDdata object
        '''

        # extract info from charmm
        psf_types   = charmm.psf.atoms[['ID', 'Type']].set_index('ID').to_dict()['Type']
        #print(psf_types)

        # substitute atoms numbers with charmm atom types
        impros     = charmm.psf.impropers.copy()
        impros['atom1'] = impros['atom1'].map(psf_types)
        impros['atom2'] = impros['atom2'].map(psf_types)
        impros['atom3'] = impros['atom3'].map(psf_types)
        impros['atom4'] = impros['atom4'].map(psf_types)
        impros['atuple'] = list(zip(impros.atom1, impros.atom2, impros.atom3))
        impros.drop(columns=['atom1', 'atom2', 'atom3', 'atom4'], inplace=True)
        #print(impros)

        # build translation dict
        btypes = impros.copy().drop_duplicates(inplace=False)
        #print(btypes)
        btypes['ID'] = np.arange(1, len(btypes)+1)
        btypes.set_index('atuple', inplace=True)
        btypeToInt = btypes.to_dict()['ID']
        #print(btypeToInt)
        btypes.reset_index(inplace=True)

        # final table
        impros['ID'] = np.arange(1, len(impros)+1)
        impros['Type'] = impros['atuple'].map(btypeToInt)
        impros.drop(columns=['atuple'], inplace=True)
        impros['Atom1'] = charmm.psf.impropers.copy()['atom1']
        impros['Atom2'] = charmm.psf.impropers.copy()['atom2']
        impros['Atom3'] = charmm.psf.impropers.copy()['atom3']
        impros['Atom4'] = charmm.psf.impropers.copy()['atom4']

        impros.rename(columns={'ID':'iID', 'Type':'iType'}, inplace=True)
        #print(impros)

        super(ImpropersDF, self).__init__(impros)
        #print(self.dtypes)


#===================================================================

class PairCoeffs(ForceField):
    def __init__(self,data=None, dtype=None, copy=False):
        super(PairCoeffs, self).__init__(data=data, columns=['aType','aType2', 'epsilon', 'sigma', 'epsilon1_4', 'sigma1_4'], dtype=dtype, copy=copy)

    def setFromNAMD(self, charmm, mass):
        ''' Extracts info from PRM and PSF objects into self.

        Parameter
        -----------------
        charmm : NAMDdata
            NAMDdata object

        mass : MassDF
            AtomsDF object associateed with these PairCoeffs
        '''

        # extract info from charmm
        psf_types   = charmm.psf.atoms[['ID', 'Type']].set_index('ID').to_dict()['Type']
        #print(psf_types)

        # substitute atoms numbers with charmm atom types
        nonbonded       = mass.copy()
        nonbonded['types'] = nonbonded.aType.map(psf_types)
        nonbonded.drop(columns=['Mass'], inplace=True)
        #print(nonbonded)

        prmFF = charmm.prm.nonbonded.getCoeffs()
        #print("Aqui el DF: ",prmFF)

        # add charge and energy to atoms
        nonbonded['epsilon'] = nonbonded.types.map(prmFF.epsilon.to_dict())
        nonbonded['sigma'] = nonbonded.types.map(prmFF.Rmin2.to_dict())
        nonbonded['epsilon1_4'] = nonbonded.types.map(prmFF.epsilon.to_dict())
        nonbonded['sigma1_4'] = nonbonded.types.map(prmFF.Rmin2.to_dict())
        nonbonded.drop(columns=['types'], inplace=True)
        
        nonbonded.rename(columns={'aID':'aType'}, inplace=True)
        #Columna nueva para el Lennard-Jones y reorganizacion columna
        nonbonded['aType2'] = nonbonded['aType']
        nonbonded = nonbonded[['aType','aType2', 'epsilon','sigma','epsilon1_4','sigma1_4']]
        print("aqui el nonbonded, la tabla:",nonbonded)

        super(PairCoeffs, self).__init__(nonbonded)
        #print("\nPairCoeffs Nans:\n",nonbonded.isna().sum())
        #print(self.dtypes)



class AngleCoeffs(ForceField):
    def __init__(self,data=None, dtype=None, copy=False):
        super(AngleCoeffs, self).__init__(data =data, columns=['anType', 'Ktheta', 'Theta0', 'Kub', 'S0'], dtype=dtype, copy=copy)

    def setFromNAMD(self, charmm, angles):
        ''' Extracts info from PRM and PSF objects into self.

        Parameter
        -----------------
        charmm : NAMDdata
            NAMDdata object

        angles : AnglesDF
            AnglesDF object associateed with these AngleCoeffs
        '''

        # extract info from charmm
        psf_types   = charmm.psf.atoms[['ID', 'Type']].set_index('ID').to_dict()['Type']
        #print(psf_types)

        # substitute atoms numbers with charmm atom types
        angles     = angles.copy()
        angles.Atom1 = angles.Atom1.map(psf_types)
        angles.Atom2 = angles.Atom2.map(psf_types)
        angles.Atom3 = angles.Atom3.map(psf_types)
        angles['atuple'] = list(zip(angles.Atom1, angles.Atom2, angles.Atom3))
        angles.drop(columns=['anID', 'Atom1', 'Atom2', 'Atom3'], inplace=True)
        angles.drop_duplicates(inplace=True)

        prmFF = charmm.prm.angles.getCoeffs()

        # add Kb and b0 to bonds
        angles['Ktheta'] = angles.atuple.map(prmFF.Ktheta.to_dict())
        #angles['Ktheta'] = angles.atuple.apply(findWithX, args=(prmFF.Ktheta.to_dict(),) )
        angles['Theta0'] = angles.atuple.map(prmFF.Theta0.to_dict())
        angles['Kub'] = angles.atuple.map(prmFF.Kub.to_dict())
        angles['S0'] = angles.atuple.map(prmFF.S0.to_dict())
        angles.drop(columns=['atuple'], inplace=True)
        #print(angles)
        #print(angles.isna().sum())

        angles.fillna(0.0, inplace=True)

        angles.rename(columns={'ID':'anType'}, inplace=True)

        super(AngleCoeffs, self).__init__(angles)
        #print("\nAngleCoeffs Nans:\n",angles.isna().sum())
        #print("Nans:\n",self)



class BondCoeffs(ForceField):
    def __init__(self,data=None, dtype=None, copy=False):
        super(BondCoeffs, self).__init__(data=data, columns=['bType','Spring_Constant','Eq_Length'], dtype=dtype, copy=copy)

    def setFromNAMD(self, charmm, bonds):
        ''' Extracts info from PRM and PSF objects into self.

        Parameter
        -----------------
        charmm : NAMDdata
            NAMDdata object

        bonds : BondsDF
            BondsDF object associateed with these BondCoeffs
        '''

        # extract info from charmm
        psf_types   = charmm.psf.atoms[['ID', 'Type']].set_index('ID').to_dict()['Type']
        #print(psf_types)

        # substitute atoms numbers with charmm atom types
        bonds     = bonds.copy()
        bonds.Atom1 = bonds.Atom1.map(psf_types)
        bonds.Atom2 = bonds.Atom2.map(psf_types)
        bonds['atuple'] = list(zip(bonds.Atom1, bonds.Atom2))
        bonds.drop(columns=['bID', 'Atom1', 'Atom2'], inplace=True)
        bonds.drop_duplicates(inplace=True)
        #print(bonds)

        prmFF = charmm.prm.bonds.getCoeffs()


        # add Kb and b0 to bonds
        bonds['Spring_Constant'] = bonds.atuple.map(prmFF.Kb.to_dict())
        bonds['Eq_Length'] = bonds.atuple.map(prmFF.b0.to_dict())
        bonds.drop(columns=['atuple'], inplace=True)
        #print(bonds)
        #print("\nBondCoeffs Nans:\n",bonds.isna().sum())

        bonds.rename(columns={'ID':'bType'}, inplace=True)

        super(BondCoeffs, self).__init__(bonds)
        #print(self.dtypes)


class DihedralCoeffs(ForceField):
    def __init__(self,data=None, dtype=None, copy=False):
        super(DihedralCoeffs, self).__init__(data=data, columns=['dType', 'Kchi', 'n', 'delta'], dtype=dtype, copy=copy)

    def setFromNAMD(self, charmm, dihedrals):
        ''' Extracts info from PRM and PSF objects into self.

        Parameter
        -----------------
        charmm : NAMDdata
            NAMDdata object

        dihedrals : DihedralsDF
            AnglesDF object associateed with these AngleCoeffs
        '''

        # extract info from charmm
        psf_types   = charmm.psf.atoms[['ID', 'Type']].set_index('ID').to_dict()['Type']

        # substitute atoms numbers with charmm atom types

        dihedrals     = dihedrals.copy()
        dihedrals.Atom1 = dihedrals.Atom1.map(psf_types)
        dihedrals.Atom2 = dihedrals.Atom2.map(psf_types)
        dihedrals.Atom3 = dihedrals.Atom3.map(psf_types)
        dihedrals.Atom4 = dihedrals.Atom4.map(psf_types)
        dihedrals['atuple'] = list(zip(dihedrals.Atom1, dihedrals.Atom2, dihedrals.Atom3, dihedrals.Atom4))
        dihedrals.drop(columns=['dID', 'Atom1', 'Atom2', 'Atom3', 'Atom4'], inplace=True)
        dihedrals.drop_duplicates(inplace=True)
        #print(dihedrals)

        prmFF = charmm.prm.dihedrals.getCoeffs()

        # add Kb and b0 to bonds
        #dihedrals['Kchi'] = dihedrals.atuple.map(prmFF.Kchi.to_dict())
        #dihedrals['n'] = dihedrals.atuple.map(prmFF.n.to_dict())
        #dihedrals['delta'] = dihedrals.atuple.map(prmFF.delta.to_dict())
        dihedrals['Kchi'] = dihedrals.atuple.apply(findWithX, args=(prmFF.Kchi.to_dict(),) )
        dihedrals['n'] = dihedrals.atuple.apply(findWithX, args=(prmFF.n.to_dict(),) )
        dihedrals['delta'] = dihedrals.atuple.apply(findWithX, args=(prmFF.delta.to_dict(),) )
        #print(dihedrals)
        dihedrals.drop(columns=['atuple'], inplace=True)
        #print("\nDihedralCoeffs Nans:\n",dihedrals.isna().sum())

        dihedrals = dihedrals.dropna()
        dihedrals['Type'] = dihedrals.index = np.arange(1, len(dihedrals)+1)
        dihedrals['Weighting_Factor'] = float(random.randint(0,2)/2)    #Is given randomly for now

        dihedrals.rename(columns={'ID':'dType'}, inplace=True)
        del dihedrals['Type']
        
        super(DihedralCoeffs, self).__init__(dihedrals.astype({'Kchi':float, 'n':int, 'delta':int,'Weighting_Factor':float}))


class ImproperCoeffs(ForceField):
    def __init__(self,data=None, dtype=None, copy=False):
        super(ImproperCoeffs, self).__init__(data=data, columns=['iType', 'Kchi', 'n', 'delta'], dtype=dtype, copy=copy)

    def setFromNAMD(self, charmm, impropers):
        ''' Extracts info from PRM and PSF objects into self.

        Parameter
        -----------------
        charmm : NAMDdata
            NAMDdata object

        impropers : ImpropersDF
            AnglesDF object associateed with these AngleCoeffs
        '''

        # extract info from charmm
        psf_types   = charmm.psf.atoms[['ID', 'Type']].set_index('ID').to_dict()['Type']
        #print(psf_types)

        # substitute atoms numbers with charmm atom types

        impropers     = impropers.copy()
        impropers.Atom1 = impropers.Atom1.map(psf_types)
        impropers.Atom2 = impropers.Atom2.map(psf_types)
        impropers.Atom3 = impropers.Atom3.map(psf_types)
        impropers.Atom4 = impropers.Atom4.map(psf_types)
        impropers['atuple'] = list(zip(impropers.Atom1, impropers.Atom2, impropers.Atom3, impropers.Atom4))
        impropers.drop(columns=['iID', 'Atom1', 'Atom2', 'Atom3', 'Atom4'], inplace=True)
        impropers.drop_duplicates(inplace=True)
        #print(impropers)

        prmFF = charmm.prm.impropers.getCoeffs()
        #print(prmFF)

        # add Kb and b0 to bonds
        impropers['Kpsi'] = impropers.atuple.apply(findWithX, args=(prmFF.Kpsi.to_dict(),) )
        impropers['psi0'] = impropers.atuple.apply(findWithX, args=(prmFF.psi0.to_dict(),) )
        #impropers['Kpsi'] = impropers.atuple.map(prmFF.Kpsi.to_dict())
        #impropers['psi0'] = impropers.atuple.map(prmFF.psi0.to_dict())
        #print(impropers)
        impropers.drop(columns=['atuple'], inplace=True)
        #print("\nImproperCoeffs Nans:\n",impropers.isna().sum())

        impropers['Type'] = impropers.index = np.arange(1, len(impropers)+1)

        impropers.rename(columns={'ID':'iType'}, inplace=True)
        del impropers['Type']
        
        super(ImproperCoeffs, self).__init__(impropers)
 

#=================================================================== 

class LammpsData():
    ''' Holds LAMMPS data.
        See "Format of a data file" section in https://lammps.sandia.gov/doc/read_data.html
    '''

    def __init__(self,file=None):

        '''
        self['Ellipsoids']  = EllipsoidsDF()
        selflines       = LinesDF()
        self.bodies      = BodiesDF()
        self.triangles   = TrianglesDF()
        self.bodies   = BodiesDF()
        '''
        
        #forceFieldSectio composi
        self.forceField = ForceFieldData()
        
        # atom-property sections
        self.atomproperty = AtomPropertyData()
        
        # molecular topology sections
        self.topologia = MolecularTopologyData()
        
        self.region = Box()
        
        if file:
            self.read(file)
        

    def read(self, filename):
    
        #Abrir el archivo para leer datos
        arch = open(filename, 'r')
        serial = 1
        ind = 0
        num = []
        tipo = []
        caja = []
        
        #l = LammpsData()
               
        #Prepare cycle to enhance process by removing CAPS and spaces to get DF keyword. 
        keywords = ["Masses","Pair Coeffs","Bond Coeffs","Angle Coeffs","Dihedral Coeffs",
                    "Improper Coeffs","Atoms","Velocities","Bonds","Angles","Dihedrals","Impropers"]

        for linea in arch:
            valid = True
            
            key = linea.strip()
            data = []
            #Si encuentro la seccion deseada en el archivo
            if key in keywords:
            
                #Comienzo a almacenar la informacion
                linea = next(arch)
                linea = next(arch)   
                       
                #Almacenar la informacion mientras exista
                while valid:      
                    data.append(linea.split())
                    #Prevenir error de iteracion cuando el archivo se acaba
                    try: linea = next(arch)
                    except(StopIteration):
                        valid = False
                        break
                    #Si se acaban los datos del parametro, vamos al proximo
                    if linea.strip() == '':
                        valid = False
                       
                #Almacenar toda la data de esa seccion
                if key == 'Masses': self.atomproperty.masses.add(data)
                if key == 'Pair Coeffs': self.forceField.pairCoeffs.add(data)
                if key == 'Bond Coeffs': self.forceField.bondCoeffs.add(data)
                if key == 'Angle Coeffs': self.forceField.angleCoeffs.add(data)                
                if key == 'Dihedral Coeffs': self.topologia.dihedralCoeffs.add(data)
                if key == 'Improper Coeffs':self.forceField.improperCoeffs.add(data)
                if key == 'Atoms': self.atomproperty.atoms.add(data)                
                if key == 'Velociteies': self.atomproperty.velocities.add(data)   
                if key == 'Bonds': self.topologia.bonds.add(data)
                if key == 'Angles': self.topologia.angles.add(data)                
                if key == 'Dihedrals': self.topologia.dihedrals.add(data)
                if key == 'Impropers':self.topologia.impropers.add(data)                
              
        arch.close()                
           

    def loadNAMDdata(self, charmm):
        ''' loads data from NAMDdata object into self.'''
        #print("loadNAMDdata=",charmm.psf.dihedrals)
        
        #AtomPropertyData
        self.atomproperty.setFromNAMD(charmm)
        
        # MolecularTopologyData
        self.topologia.setFromNAMD(charmm)
       
        # ForceFieldData
      
        self.forceField.setFromNAMD(charmm, self.atomproperty,self.topologia)
        
        # MolecularTopologyData
        self.region.setFromNAMD(charmm)
       
        '''
        self.pairCoeffs.setFromNAMD(charmm, self.masses)
        self.bondCoeffs.setFromNAMD(charmm, self.bonds)
        self.angleCoeffs.setFromNAMD(charmm, self.angles)
        self.dihedralCoeffs.setFromNAMD(charmm, self.dihedrals)
        self.improperCoeffs.setFromNAMD(charmm, self.impropers)
        '''

    def loadWolffia(self, wolffia):
        '''
        Converts a WolffiaState to a LammpsData self object.

        @param: mix a Wolffia Mixture object.
        '''
        from granules.structure.NAMDdata import NAMDdata

        charmm = NAMDdata()
        charmm.loadWolffia(wolffia)
        self.loadNAMDdata(charmm)
 
        return self


    def writeConf(self, filename):
        import sys

        cfile = open(filename, "w")


        # Sección automática
        overview_section = '{0:>12}  atoms\n{1:>12}  bonds\n{2:>12}  angles\n{3:>12}  dihedrals\n{4:>12}  impropers\n\n{5:>12}  atom types\n' \
                           '{6:>12}  bond types\n{7:>12}  angle types\n{8:>12}  dihedral types\n{9:>12}  improper types\n\n' \
                           .format( \
                                len(self.atomproperty.atoms), len(self.topologia.bonds), len(self.topologia.angles), len(self.topologia.dihedrals), len(self.topologia.impropers), len(self.atomproperty.masses), \
                                len(self.forceField.bondCoeffs), len(self.forceField.angleCoeffs), len(self.forceField.dihedralCoeffs), len(self.forceField.improperCoeffs)) 

        try:
            maxminstr = [str(x) for x in self.region.maxsMins]
            box_section =  ' '.join(maxminstr[:2])  + ' xlo xhi\n' + \
                           ' '.join(maxminstr[2:4]) + ' ylo yhi\n' + \
                           ' '.join(maxminstr[4:])  + ' zlo zhi\n'
        except:
            # Improve and ask from where we get this data
            box_section =  ' ' + str(self.atomproperty.atoms['x'].min()-2) + ' ' + str(self.atomproperty.atoms['x'].max()+2) + ' xlo xhi\n' + \
                           ' ' + str(self.atomproperty.atoms['y'].min()-2) + ' ' + str(self.atomproperty.atoms['y'].max()+2) + ' ylo yhi\n' + \
                           ' ' + str(self.atomproperty.atoms['z'].min()-2) + ' ' + str(self.atomproperty.atoms['z'].max()+2) + ' zlo zhi\n'   
                       
        # Header
        cfile.write("LAMMPS Description\n\n")
        cfile.write(overview_section + box_section)


        #Masses
        if len(self.atomproperty.masses) > 0:
            cfile.write('\nMasses\n\n')
            cfile.write(self.atomproperty.masses.to_string(index=False, columns=self.atomproperty.masses.columns, header=False))
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No atom mass values to write. Simulation unlikely to run with this file.")


        #Pair Coeffs
        if len(self.forceField.pairCoeffs) > 0:
            cfile.write('\nPair Coeffs\n\n')
            #print("Pair Coeffs:", self.pairCoeffs.columns)
            #[cfile.write('{:>3d}{:>12}{:>12}{:>12}{:>12}\n'.format(row['ID'], row['Charge'], row['Energy'],
            # row['Charge'], row['Energy'])) for index, row in self['Pair Coeffs'].iterrows()]
            cfile.write(self.forceField.pairCoeffs.to_string(index=False, columns=self.forceField.pairCoeffs.columns, header=False))        
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No Pair coefficients to write.\n")


        #Bond Coeffs
        if len(self.forceField.bondCoeffs) > 0:
            cfile.write('\nBond Coeffs\n\n')
            cfile.write(self.forceField.bondCoeffs.to_string(index=False, columns=self.forceField.bondCoeffs.columns, header=False))        
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No bond coefficients to write.\n")


        #Angle Coeffs
        if len(self.forceField.angleCoeffs) > 0:
            cfile.write('\nAngle Coeffs\n\n')
            cfile.write(self.forceField.angleCoeffs.to_string(index=False, columns=self.forceField.angleCoeffs.columns, header=False))        
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No angle coefficients to write.\n")


        #Dihedral Coeffs
        if len(self.forceField.dihedralCoeffs) > 0:
            cfile.write('\nDihedral Coeffs\n\n')
            cfile.write(self.forceField.dihedralCoeffs.to_string(index=False, columns=self.forceField.dihedralCoeffs.columns, header=False))        
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No dihedral coefficients to write.\n")


        #Improper Coeffs
        if len(self.forceField.improperCoeffs) > 0:
            cfile.write('\nImproper Coeffs\n\n') 
            cfile.write(self.forceField.improperCoeffs.to_string(index=False, columns=self.forceField.improperCoeffs.columns, header=False))        
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No improper coefficients to write.\n")


        #Atoms
        cfile.write('\nAtoms\n\n') 
        #cfile.write(self.atomproperty.atoms.to_string(index=False, columns=self.atomproperty.atoms.columns, header=False))
        cfile.write(self.atomproperty.atoms.to_string(index=False, columns=self.atomproperty.atoms.columns, header=False))
        cfile.write("\n")


        #Velocities
        if len(self.atomproperty.velocities) > 0:
            cfile.write('\nVelocities\n\n')
            cfile.write(self.atomproperty.velocities.to_string(index=False, columns=self.atomproperty.velocities.columns, header=False))        
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No velocities to write.\n")


        #Bonds
        if len(self.topologia.bonds) > 0:
            cfile.write('\nBonds\n\n')
            cfile.write(self.topologia.bonds.to_string(index=False, columns=self.topologia.bonds.columns, header=False))
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No bonds to write.\n")

        
        #Angles
        if len(self.topologia.angles) > 0:
            cfile.write('\nAngles\n\n') 
            cfile.write(self.topologia.angles.to_string(index=False, columns=self.topologia.angles.columns, header=False))
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No angles to write.\n")


        #Dihedrals
        if len(self.topologia.dihedrals) > 0:
            cfile.write('\nDihedrals\n\n') 
            cfile.write(self.topologia.dihedrals.to_string(index=False, columns=self.topologia.dihedrals.columns, header=False))
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No dihedrals to write.\n")


        #Impropers
        if len(self.topologia.impropers) > 0:
            cfile.write('\nImpropers\n\n') 
            cfile.write(self.topologia.impropers.to_string(index=False, columns=self.topologia.impropers.columns, header=False))
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No impropers to write.\n")
    
    
    def charmmForce(self):
        '''Hace una llamada a la funcion charmmForce() de la clase forceField() 
            para poder printiar los datos.'''
        
        print("LammpsData.charmmForce()")
        self.forceField.charmmForce(self.atomproperty,self.topologia)

    def append(self,other):
        '''Une dos objetos de LammpsData, sus dataframes individuales'''
        # OH YEAH
        self.forceField.angleCoeffs = self.forceField.angleCoeffs.append(other.forceField.angleCoeffs)
        self.forceField.bondCoeffs = self.forceField.bondCoeffs.append(other.forceField.bondCoeffs)
        self.forceField.dihedralCoeffs = self.forceField.dihedralCoeffs.append(other.forceField.dihedralCoeffs)
        self.forceField.improperCoeffs = self.forceField.improperCoeffs.append(other.forceField.improperCoeffs)
        self.forceField.pairCoeffs = self.forceField.pairCoeffs.append(other.forceField.pairCoeffs)
        #Oh YEaSH
        self.atomproperty.atoms = self.atomproperty.atoms.append(other.atomproperty.atoms)
        self.atomproperty.velocities = self.atomproperty.velocities.append(other.atomproperty.velocities)
        self.atomproperty.masses = self.atomproperty.masses.append(other.atomproperty.masses)
        #yeyeyeye
        self.topologia.angles = self.topologia.angles.append(other.topologia.angles)
        self.topologia.bonds = self.topologia.bonds.append(other.topologia.bonds)
        self.topologia.dihedrals = self.topologia.dihedrals.append(other.topologia.dihedrals)
        self.topologia.impropers = self.topologia.impropers.append(other.topologia.impropers)

        
    def copy(self):#modifica
        ld = LammpsData()
        # OH YEAH
        ld.forceField.angleCoeffs = self.forceField.angleCoeffs.copy()
        ld.forceField.bondCoeffs = self.forceField.bondCoeffs.copy()
        ld.forceField.dihedralCoeffs = self.forceField.dihedralCoeffs.copy()
        ld.forceField.improperCoeffs = self.forceField.improperCoeffs.copy()
        ld.forceField.pairCoeffs = self.forceField.pairCoeffs.copy()
        # OH YEAHHH
        ld.atomproperty.atoms = self.atomproperty.atoms.copy()
        ld.atomproperty.velocities = self.atomproperty.velocities.copy()
        ld.atomproperty.masses = self.atomproperty.masses.copy()
        # OH YEAHHH
        ld.topologia.angles = self.topologia.angles.copy()
        ld.topologia.bonds = self.topologia.bonds.copy()
        ld.topologia.dihedrals = self.topologia.dihedrals.copy()
        ld.topologia.impropers = self.topologia.impropers.copy()
        
        return ld
    
    def selectAtom(self,atomNumber):
        '''Funcion que elimina un tipo de atomo deseado del dataframe'''
        
        # AtompropertyData 
        #atoms
        self.atomproperty.atoms.drop(self.atomproperty.atoms[self.atomproperty.atoms.aType == atomNumber].index, inplace=True)#elimina el atomo indicado
        self.atomproperty.atoms.reset_index(inplace=True)#reset el index, crea una copia del index
        self.atomproperty.atoms.drop(['index','aID'],axis = 1, inplace=True)#elimina la copia y la columna afectada
        self.atomproperty.atoms.insert(0, "aID", [*range(1,len(self.atomproperty.atoms)+1)], True)#inserta columna nueva en la afectada("aID")
        #velocity
        self.atomproperty.velocities.drop(self.atomproperty.velocities.tail(166-len(self.atomproperty.atoms)).index,inplace = True) #elimina velocidades dependiendo a los atomos eliminados
        #masses
         #chequea
         
        #TopologiaData
        #bonds
        #Cambio para que busque en los Atomos enves de Type
        self.topologia.bonds = self.topologia.bonds.ix[self.topologia.bonds['Atom1'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topologia.bonds = self.topologia.bonds.ix[self.topologia.bonds['Atom2'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topologia.bonds  = self.topologia.bonds.reset_index(drop=True)#resetea el indice y elimina la copia
        self.topologia.bonds['bID'] = [i for i in range(1,len(self.topologia.bonds)+1 )]#resetea el indice en la columna 'bType'
        #angles
        self.topologia.angles = self.topologia.angles.ix[self.topologia.angles['Atom3'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topologia.bonds = self.topologia.bonds.ix[self.topologia.bonds['Atom1'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topologia.bonds = self.topologia.bonds.ix[self.topologia.bonds['Atom2'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topologia.angles  = self.topologia.angles.reset_index(drop=True)#resetea el indice y elimina la copia
        self.topologia.angles['anID'] = [i for i in range(1,len(self.topologia.angles)+1 )]#resetea el indice en la columna 'anType'
        #improper######modifica desde aqui
        self.topologia.impropers = self.topologia.impropers.ix[self.topologia.impropers['iType'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topologia.impropers  = self.topologia.impropers.reset_index(drop=True)#resetea el indice y elimina la copia
        self.topologia.impropers['iID'] = [i for i in range(1,len(self.topologia.impropers)+1 )]#resetea el indice en la columna 'iType'
        #dihedral
        self.topologia.dihedrals = self.topologia.dihedrals.ix[self.topologia.dihedrals['dType'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topologia.dihedrals  = self.topologia.dihedrals.reset_index(drop=True)#resetea el indice y elimina la copia
        self.topologia.dihedrals['dID'] = [i for i in range(1,len(self.topologia.dihedrals)+1 )]#resetea el indice en la columna 'dType'
        
        #ForceFieldData 
        #self.forceField.angleCoeffs
            #Falta
        
        
      
       

        
        '''
        #forceFieldSectio composi
        self.forceField = ForceFieldData()
        # molecular topology sections
        self.topologia = MolecularTopologyData()
        '''
      
_AVOGRADRO_CONSTANT_ = 6.02214129e+23

class Region:
    '''
        Container for the Mixture.
    '''
    
    DENSITY = {None:[0,0,1],"WATER": [0.9970479,298.15,18.01528], "Chloroform":[1.483,298.15,119.38], "Chloroform (4 atoms)": [1.483,298.15,119.38], "Acetone": [0.786,293.15,58.08], "THF": [0.886,293.15,72.11], "Argon": [5.537, 0.03049,39.948], "Methanol": [0.791,293.15,32.04], "Selected Molecule": [0.9970479,298.15,18.01528]}
    DI = 0
    TI = 1
    MI = 2

    def __init__(self, cid=None):
        '''
            id (int): identifier (LAMPS style)
        '''
        # Lammps-style variables
        self.id = cid
        
    def lammpsCommand(self, mixture):
        ''' to be done '''
        pass

    def setId(self, cid): self.id = cid
    
    def amountSolventMolecules(self, solvent):
            #import numpy
            global _AVOGRADRO_CONSTANT_
            D = self.DENSITY[solvent][self.DI]
            MM = self.DENSITY[solvent][self.MI]
            V = self.volume() / 1E24
            if V == 0: return 0
            
            #print "computeMolecules ", D,V,MM,D * V / MM * _AVOGRADRO_CONSTANT_
            n = int(D * V / MM * _AVOGRADRO_CONSTANT_)
    
            return n


# ===========================================
class Box(Region):
    '''
        Defines a box with rectangular anlges (similar to block in LAMMPS)
    '''
    def __init__(self, cid=None, maxsMins=None):
        '''        
            maxsMins = (xmin, xmax, ymin, ymax, zmin, zmax) or None
        '''
        super(Box, self).setId(cid)
    
        self.setMinsMaxs(maxsMins)
    
    def getMaxsMins(self): return self.maxsMins
    
    def getCenter(self):
        cellBasisVectors = self.getCellBasisVectors()
        x = (cellBasisVectors[0][0]+cellBasisVectors[1][0]+cellBasisVectors[2][0])/2
        y = (cellBasisVectors[0][1]+cellBasisVectors[1][1]+cellBasisVectors[2][1])/2
        z = (cellBasisVectors[0][2]+cellBasisVectors[1][2]+cellBasisVectors[2][2])/2
        return [x,y,z]

    def getFaces(self):
        cellBasisVectors = self.getCellBasisVectors()
        left = self.cellOrigin[0]
        right = left + cellBasisVectors[0][0]
        bottom = self.cellOrigin[1]
        top = bottom + cellBasisVectors[1][1]
        zNear = self.cellOrigin[2]
        zFar = zNear + cellBasisVectors[2][2]
            
        return [left, right, bottom, top, zNear, zFar]

    def lammpsCommand(self):
        return "region {:d} style block {}".format(self.id, ' '.join(self.maxsMins))
    
    def setFromNAMD(self, charmm):
        ''' Extracts info from NAMD.data.PBC object thet is assumed to represent a box

        Parameter
        -----------------
        charmm : NAMDdata
            NAMDdata object
        '''
        
        try:
            self.setMinsMaxs([charmm.pbc.cellOrigin[0],
                              charmm.pbc.cellBasisVector1[0] + charmm.pbc.cellOrigin[0],
                              charmm.pbc.cellOrigin[1],
                              charmm.pbc.cellBasisVector2[1] + charmm.pbc.cellOrigin[1],
                              charmm.pbc.cellOrigin[2],
                              charmm.pbc.cellBasisVector3[2] + charmm.pbc.cellOrigin[2]
                              ])
        except: 
            self.setMinsMaxs(None)

    def loadFromDump(self, filename):
        ''' Extracts info from LAMMPS dump file that is assumed to represent a box

        Parameter
        -----------------
        filename : LAMMPS dump file
        '''
        dump = open(filename, "r")
        for linea in dump:
            if linea[:26] == "ITEM: BOX BOUNDS pp pp pp":
                mismaxsstr = linea.readline().split(' ').append(
                        linea.readline().split(' ')).append(
                        linea.readline().split(' '))
                self.setMinsMaxs([float(x) for x in mismaxsstr])
                break
        dump.close()

    def setMinsMaxs(self, maxsMins):
        self.maxsMins = maxsMins
        
    def volume(self):
            '''
            volume
            '''
            return (self.maxsMins[1]-self.maxsMins[0]) * \
                (self.maxsMins[3]-self.maxsMins[2]) * \
                (self.maxsMins[5]-self.maxsMins[4])



if __name__ == "__main__":  # tests
    from NAMDdata import NAMDdata
    
    l = LammpsData()
    l.read('data.peptide')
    
    #Imprimir todos los datos    
    #input(l.__dict__)

    '''
    ch = NAMDdata()
    ch.readFiles("2rvd_autopsf.pdb", "2rvd_autopsf.psf", "par_all36_prot.prm")

    l.loadNAMDdata(ch)
    l.writeConf("2rvd.data")
    '''
