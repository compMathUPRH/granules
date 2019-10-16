# -*- coding: utf-8 -*-
#---------------------------------------------------------------------------
"""
  LAMMPSdata.py
  First version, October, 2011


  Part of granules Version 0.1.0, October, 2019
    Copyright 2019: José O.  Sotero Esteva, Lyxaira M. Glass Rivera, 
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
        #print("\nVariableColumnsSections add ", type(self))
        #Conocer cuantas columnas deben haber
        cont=len(data[0])
        columns = list(self.columns)
        col_num=len(columns)

        #Guardar los nombres adecuados de las columnas
        while cont != col_num:
            col = 'Coeff' + str(col_num)
            columns.append(col)
            col_num += 1
       
        #Separar la data para almacenar en diccionario
        dtable = {col:[] for col in columns}
        #print("VariableColumnsSections columns ", columns)
        #print("VariableColumnsSections data ", data[0])
          
        for d in data:
            for col in columns:
                try:
                    dtable[col].append(self[col].dtype.type(d.pop(0)))
                except KeyError:
                    dtable[col].append(float(d.pop(0)))


        #Crear dataframe con toda la informacion
        fr = pd.DataFrame(dtable, copy=False, columns=columns)
       
        #Ajustar indices para que no comiencen en 0
        fr.index = np.arange(1, len(fr) + 1)
       
        #Reescribir DataFrame original con creado
        if self.empty:
            super().__init__(data=fr)
        else:
            super().__init__(data=self.append(fr, ignore_index=True))
        #print("VariableColumnsSections types ", self.dtypes)


#===================================================================

class AtomPropertySection(LammpsBodySection):
    pass

class MolecularTopologySections(LammpsBodySection):
    pass

class ForceFieldSection(LammpsBodySection):
    pass

#===================================================================
   
class AtomsDF(AtomPropertySection):
    def __init__(self,data=None, dtype=None, copy=False):
        dtypes = {'aID':[0], 'Mol_ID':[0], 'Type':[0], 'Q':[0.0], 
                  'X':[0.0], 'Y':[0.0], 'Z':[0.0], 'Nx':[0], 'Ny':[0], 'Nz':[0]}   
        super(AtomsDF, self).__init__(data=dtypes, copy=copy, columns=['aID', 'Mol_ID', 'Type', 'Q', 
                  'X', 'Y', 'Z', 'Nx', 'Ny', 'Nz'])
        super(AtomsDF, self).__init__(self.drop([0]))

    def setFromPSF(self, charmm):
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
        sel_pdb     = charmm.pdb[['ID','x','y','z']].set_index('ID')
        sel         = sel_pdb.join(sel_psf)
        sel['Type'] = sel_psf['Type'].map(charmmTypeToInt)
        sel.reset_index(inplace=True)
        sel         .rename(columns={"Charge":"Q",'ID':'aID'}, inplace=True)

        # add remining columns
        sel['Mol_ID'] = np.ones((len(sel), 1), dtype=np.int8)
        sel['Nx']     = np.zeros((len(sel), 1))
        sel['Ny']     = np.zeros((len(sel), 1))
        sel['Nz']     = np.zeros((len(sel), 1))
        sel['ID']     = np.arange(1, len(sel)+1)

        # rearrange columns
        sel = sel[['aID', 'Mol_ID', 'Type', 'Q', 'x', 'y', 'z', 'Nx', 'Ny', 'Nz']]

        #sel.reset_index(inplace=True)
        #print("sel = ", sel.dtypes)
        super(AtomsDF, self).__init__(sel.astype({
                     'Mol_ID' :int,
                     'Type' :int,
                     'Q' :float,
                     'x' :float,
                     'y' :float,
                     'z' :float,
                     'Nx' :int,
                     'Ny' :int,
                     'Nz' : int
                    }))
        #print(self.dtypes)


class MassesDF(AtomPropertySection):
    def __init__(self,data=None, dtype=None, copy=False):
        if data  is None:
            dtypes = {'aID':[0], 'Mass':[0.0]}
            super(MassesDF, self).__init__(data=dtypes, copy=copy, columns=dtypes.keys())
            super(MassesDF, self).__init__(self.drop([0]))
 
    def setFromPSF(self, psf_atoms, atoms):
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
        sel_self = atoms[['aID', 'Type']].set_index('aID').copy()
        sel      = sel_self.join(sel_psf).drop_duplicates().reset_index()

        # rename columns
        sel      .drop(columns='Type', inplace=True)
        #sel      .rename(columns={"Type":"mID"}, inplace=True)
        #print(sel.dtypes)

        super(MassesDF, self).__init__(sel)


class VelocitiesDF(AtomPropertySection):
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

class AnglesDF(MolecularTopologySections):
    def __init__(self,data=None, dtype=None, copy=False):
        dtypes = {'anID':[0], 'anType':[0], 'Atom1':[0], 'Atom2':[0], 'Atom3':[0]}
        super(AnglesDF, self).__init__(data=dtypes, copy=copy, columns=dtypes.keys())
        super(AnglesDF, self).__init__(self.drop([0]))

    def setFromPSF(self, charmm):
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

class BondsDF(MolecularTopologySections):
    def __init__(self,data=None, dtype=None, copy=False):
        dtypes = {'bID':[0], 'bType':[0], 'Atom1':[0], 'Atom2':[0]}     
        super(BondsDF, self).__init__(data=dtypes, copy=copy, columns=dtypes.keys())
        super(BondsDF, self).__init__(self.drop([0]))

    def setFromPSF(self, charmm):
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
        bonds     = charmm.psf.bonds.copy()
        bonds['atom1'] = bonds['atom1'].map(psf_types)
        bonds['atom2'] = bonds['atom2'].map(psf_types)
        bonds['atuple'] = list(zip(bonds.atom1, bonds.atom2))
        bonds.drop(columns=['atom1', 'atom2'], inplace=True)
        #print(bonds)

        # build translation dict
        btypes = bonds.copy().drop_duplicates(inplace=False)
        #print(btypes)
        btypes['ID'] = np.arange(1, len(btypes)+1)
        btypes.set_index('atuple', inplace=True)
        btypeToInt = btypes.to_dict()['ID']
        #print(btypeToInt)
        btypes.reset_index(inplace=True)

        # final table
        bonds['ID'] = np.arange(1, len(bonds)+1)
        bonds['Type'] = bonds['atuple'].map(btypeToInt)
        bonds.drop(columns=['atuple'], inplace=True)
        bonds['Atom1'] = charmm.psf.bonds.copy()['atom1']
        bonds['Atom2'] = charmm.psf.bonds.copy()['atom2']

        bonds.rename(columns={'ID':'bID', 'Type':'bType'}, inplace=True)
        #print(bonds)

        super(BondsDF, self).__init__(bonds)
        #print(self.dtypes)


       
class DihedralsDF(MolecularTopologySections):
    def __init__(self,data=None, dtype=None, copy=False):
        dtypes = {'dID':[0], 'dType':[0], 'atom1':[0], 'atom2':[0], 'atom3':[0], 'atom4':[0]}
        super(DihedralsDF, self).__init__(data=dtypes, copy=copy, columns=dtypes.keys())
        super(DihedralsDF, self).__init__(self.drop([0]))

    def setFromPSF(self, charmm):
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


class ImpropersDF(MolecularTopologySections):
    def __init__(self,data=None, dtype=None, copy=False):
        dtypes = {'iID':[0], 'iType':[0], 'atom1':[0], 'atom2':[0], 'atom3':[0], 'atom4':[0]}
        super(ImpropersDF, self).__init__(data=dtypes, copy=copy, columns=dtypes.keys())
        super(ImpropersDF, self).__init__(self.drop([0]))

    def setFromPSF(self, charmm):
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

class PairCoeffs(ForceFieldSection):
    def __init__(self,data=None, dtype=None, copy=False):
        super(PairCoeffs, self).__init__(data=data, columns=['aType', 'epsilon', 'sigma', 'epsilon1_4', 'sigma1_4'], dtype=dtype, copy=copy)

    def setFromPSF(self, charmm, mass):
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
        nonbonded['types'] = nonbonded.aID.map(psf_types)
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
        #print(nonbonded)

        nonbonded.rename(columns={'ID':'aType'}, inplace=True)

        super(PairCoeffs, self).__init__(nonbonded)
        #print("\nPairCoeffs Nans:\n",nonbonded.isna().sum())
        #print(self)



class AngleCoeffs(ForceFieldSection):
    def __init__(self,data=None, dtype=None, copy=False):
        super(AngleCoeffs, self).__init__(data=data, columns=['anType', 'Ktheta', 'Theta0', 'Kub', 'S0'], dtype=dtype, copy=copy)

    def setFromPSF(self, charmm, angles):
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



class BondCoeffs(ForceFieldSection):
    def __init__(self,data=None, dtype=None, copy=False):
        super(BondCoeffs, self).__init__(data=data, columns=['bType','Spring_Constant','Eq_Length'], dtype=dtype, copy=copy)

    def setFromPSF(self, charmm, bonds):
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


class DihedralCoeffs(ForceFieldSection):
    def __init__(self,data=None, dtype=None, copy=False):
        super(DihedralCoeffs, self).__init__(data=data, columns=['dType', 'Kchi', 'n', 'delta'], dtype=dtype, copy=copy)

    def setFromPSF(self, charmm, dihedrals):
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
        #print(psf_types)

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
        #print(prmFF)

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

        super(DihedralCoeffs, self).__init__(dihedrals.astype({'Kchi':float, 'n':int, 'delta':int,'Weighting_Factor':float}))
        #print("DihedralCoeffs\n",self.dtypes)


class ImproperCoeffs(ForceFieldSection):
    def __init__(self,data=None, dtype=None, copy=False):
        super(ImproperCoeffs, self).__init__(data=data, columns=['iType', 'Kchi', 'n', 'delta'], dtype=dtype, copy=copy)

    def setFromPSF(self, charmm, impropers):
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

        super(ImproperCoeffs, self).__init__(impropers)
        #print(self)

#===================================================================
#CLASES TEMPORERAS
        
class NonbDF(MolecularTopologySections):
    def __init__(self,data=None, dtype=None, copy=False):
        super(NonbDF, self).__init__(data=data, copy=copy, columns=['ID','Type1','Force','Charge','Energy'])

class BondDF(MolecularTopologySections):
    def __init__(self,data=None, dtype=None, copy=False):
        super(BondDF, self).__init__(data=data, copy=copy, columns=['ID','Type1','Type2','Spring_Constant','Eq_Length'])

class AnglDF(MolecularTopologySections):
    def __init__(self,data=None, dtype=None, copy=False):
        super(AnglDF, self).__init__(data=data, copy=copy, columns=['ID','Type1','Type2','Type3','Spring_Constant','Eq_Angle'])

class ImprDF(MolecularTopologySections):
    def __init__(self,data=None, dtype=None, copy=False):
        super(ImprDF, self).__init__(data=data, copy=copy, columns=['ID','Atom','Types','Kpsi','Psi0'])

class DiheDF(MolecularTopologySections):
    def __init__(self,data=None, dtype=None, copy=False):
        super(DiheDF, self).__init__(data=data, copy=copy, columns=['ID','Type1','Type2','Type3','Type4','Spring_Constant','Multiplicity','Eq_Angle'])

class PDB_File(MolecularTopologySections):
    def __init__(self,data=None, dtype=None, copy=False):
        super(PDB_File, self).__init__(data=data, copy=copy, columns=['RecName','ID','Name','AltLoc','ResName','ChainID','ResSeq','iCode','x','y','z','Occupancy','TempFactor','Element','Charge'])

class PSF_File(MolecularTopologySections):
    def __init__(self,data=None, dtype=None, copy=False):
        super(PSF_File, self).__init__(data=data, copy=copy, columns=['ID','RecName','ChainID', 'ResName', 'Name', 'Type', 'Charge', 'Mass', 'Unused'])

#=================================================================== 

class LammpsData:
    ''' Holds LAMMPS data.
        See "Format of a data file" section in https://lammps.sandia.gov/doc/read_data.html
    '''

    def __init__(self):

        # atom-property sections
        self.atoms       = AtomsDF()
        self.velocities  = VelocitiesDF()
        self.masses      = MassesDF()

        '''
        self['Ellipsoids']  = EllipsoidsDF()
        selflines       = LinesDF()
        self.bodies      = BodiesDF()
        self.triangles   = TrianglesDF()
        self.bodies   = BodiesDF()
        '''

        # molecular topology sections
        self.bonds       = BondsDF()
        self.angles      = AnglesDF()
        self.dihedrals   = DihedralsDF()
        self.impropers   = ImpropersDF()

        #  force field sections
        #self.pairIJCoeffs           = PairDF()
        self.bondCoeffs             = BondCoeffs()
        self.pairCoeffs             = PairCoeffs()
        self.angleCoeffs            = AngleCoeffs()
        self.dihedralCoeffs         = DihedralCoeffs()
        self.improperCoeffs         = ImproperCoeffs()

        # class 2 force field sections
        self.bondCoeffs         = BondCoeffs()
        self.bondAngleCoeffs        = ForceFieldSection()
        self.middleBondTorsion      = ForceFieldSection()
        self.angleAngleCoeffs       = ForceFieldSection()
        self.bondBond13Coeffs       = ForceFieldSection()
        self.angleTorsionCoeffs     = ForceFieldSection()
        self.endBondTorsionCoeffs   = ForceFieldSection()
        self.angleAngleTorsionCoeffs= ForceFieldSection()
        

    def charmmBondEnergy(self):
        ''' Computes CHARMM bond energy.
            Formula: sum K * (bij - b0)**2
        '''
        bi = self.bonds.set_index('Atom1').join(
                self.atoms.set_index('aID')
             )[['bID','x', 'y', 'z']].set_index('bID')

        bj = self.bonds.set_index('Atom2').join(
                self.atoms.set_index('aID')
             )[['bID','x', 'y', 'z']].set_index('bID')

        bij = np.sqrt(np.sum((bi - bj) ** 2, axis='columns'))  # bonds lengths

        coeffs = self.bonds[['bID','bType']].set_index('bType').join(
                        self.bondCoeffs.set_index('bType')
                 ).reset_index(drop=True)
        K = coeffs[['bID','Spring_Constant']].set_index('bID').Spring_Constant
        b0 = coeffs[['bID','Eq_Length']].set_index('bID').Eq_Length

        return np.sum(K * (bij-b0)**2)

    def charmmAngleEnergy(self):
        ''' Computes CHARMM angle energy.
            Formula: sum K * (aij - a0)**2
        '''
        bi = self.angles.set_index('Atom1').join(
                self.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        bj = self.angles.set_index('Atom2').join(
                self.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        bk = self.angles.set_index('Atom3').join(
                self.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        # compute angles
        l1 = bi - bj
        l2 = bk - bj

        norm1 = np.sqrt(np.square(l1).sum(axis=1))
        norm2 = np.sqrt(np.square(l1).sum(axis=1))
        dot = l1.x * l2.x + l1.y * l2.y + l1.z * l2.z
        angles = np.arccos(dot / (norm1 * norm2))

        coeffs = self.angles[['anID','anType']].set_index('anType').join(
                        self.angleCoeffs.set_index('anType')
                 ).reset_index(drop=True)

        K = coeffs[['anID','Ktheta']].set_index('anID').Ktheta
        a0 = np.radians(coeffs[['anID','Theta0']].set_index('anID').Theta0)

        return np.sum(K * (angles-a0)**2)

    def charmmDihedralsEnergy(self):
        ''' Computes CHARMM angle energy.
            Formula: sum K * (1 + cos(n * x - d))
        '''
        bi = self.dihedrals.set_index('Atom1').join(
                self.atoms.set_index('aID')
             )[['dID','x', 'y', 'z']].set_index('dID')

        bj = self.dihedrals.set_index('Atom2').join(
                self.atoms.set_index('aID')
             )[['dID','x', 'y', 'z']].set_index('dID')

        bk = self.dihedrals.set_index('Atom3').join(
                self.atoms.set_index('aID')
             )[['dID','x', 'y', 'z']].set_index('dID')

        bl = self.dihedrals.set_index('Atom4').join(
                self.atoms.set_index('aID')
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
        print(self.dihedrals)
        coeffs = self.dihedrals[['dID','dType']].set_index('dType').join(
                        self.dihedralCoeffs.set_index('dType')
                 ).reset_index(drop=True)
        K = coeffs[['dID','Kchi']].set_index('dID').Kchi
        n = coeffs[['dID','n']].set_index('dID').n
        d = coeffs[['dID','delta']].set_index('dID').delta



        return np.sum(K * (1 + np.cos(n * angles - d)))
        '''
        '''

    def loadNAMDdata(self, charmm):
        ''' loads data from NAMDdata object into self.'''
        #print("loadNAMDdata=",charmm.psf.dihedrals)
        self.atoms.setFromPSF(charmm)
        self.velocities.setToZero(self.atoms)
        self.masses.setFromPSF(charmm.psf.atoms, self.atoms)

        # molecular topology sections
        self.bonds.setFromPSF(charmm)
        self.angles.setFromPSF(charmm)
        self.dihedrals.setFromPSF(charmm)
        self.impropers.setFromPSF(charmm)

        #  force field sections
        self.pairCoeffs.setFromPSF(charmm, self.masses)
        self.bondCoeffs.setFromPSF(charmm, self.bonds)
        self.angleCoeffs.setFromPSF(charmm, self.angles)
        self.dihedralCoeffs.setFromPSF(charmm, self.dihedrals)
        self.improperCoeffs.setFromPSF(charmm, self.impropers)

    def writeConf(self, filename):
        cfile = open(filename, "w")


        # Sección automática
        overview_section = '{0:>12}  atoms\n{1:>12}  bonds\n{2:>12}  angles\n{3:>12}  dihedrals\n{4:>12}  impropers\n\n{5:>12}  atom types\n' \
                           '{6:>12}  bond types\n{7:>12}  angle types\n{8:>12}  dihedral types\n{9:>12}  improper types\n\n' \
                           .format( \
                                len(self.atoms), len(self.bonds), len(self.angles), len(self.dihedrals), len(self.impropers), len(self.masses), \
                                len(self.bondCoeffs), len(self.angleCoeffs), len(self.dihedralCoeffs), len(self.improperCoeffs)) 

        # Improve and ask from where we get this data
        box_section =  ' ' + str(36.840194) + ' ' + str(64.211568) + ' xlo xhi\n' + \
                       ' ' + str(41.013691) + ' ' + str(68.385058) + ' ylo yhi\n' + \
                       ' ' + str(29.768095) + ' ' + str(57.139462) + ' zlo zhi\n'   \

        # Header
        cfile.write("LAMMPS Description\n\n")
        cfile.write(overview_section + box_section)


        #Masses
        cfile.write('\nMasses\n\n')
        cfile.write(self.masses.to_string(index=False, columns=self.masses.columns, header=False))
        cfile.write("\n")


        #Pair Coeffs
        cfile.write('\nPair Coeffs\n\n')
        #print("Pair Coeffs:", self.pairCoeffs.columns)
        #[cfile.write('{:>3d}{:>12}{:>12}{:>12}{:>12}\n'.format(row['ID'], row['Charge'], row['Energy'],
        # row['Charge'], row['Energy'])) for index, row in self['Pair Coeffs'].iterrows()]
        cfile.write(self.pairCoeffs.to_string(index=False, columns=self.pairCoeffs.columns, header=False))        
        cfile.write("\n")


        #Bond Coeffs
        cfile.write('\nBond Coeffs\n\n')
        cfile.write(self.bondCoeffs.to_string(index=False, columns=self.bondCoeffs.columns, header=False))        
        cfile.write("\n")


        #Angle Coeffs
        cfile.write('\nAngle Coeffs\n\n')
        cfile.write(self.angleCoeffs.to_string(index=False, columns=self.angleCoeffs.columns, header=False))        
        cfile.write("\n")


        #Dihedral Coeffs
        cfile.write('\nDihedral Coeffs\n\n')
        cfile.write(self.dihedralCoeffs.to_string(index=False, columns=self.dihedralCoeffs.columns, header=False))        
        cfile.write("\n")


        #Improper Coeffs
        cfile.write('\nImproper Coeffs\n\n') 
        cfile.write(self.improperCoeffs.to_string(index=False, columns=self.improperCoeffs.columns, header=False))        
        cfile.write("\n")


        #Atoms
        cfile.write('\nAtoms\n\n') 
        #cfile.write(self.atoms.to_string(index=False, columns=self.atoms.columns, header=False))
        cfile.write(self.atoms.to_string(index=False, columns=self.atoms.columns, header=False))
        cfile.write("\n")


        #Velocities
        cfile.write('\nVelocities\n\n')
        cfile.write(self.velocities.to_string(index=False, columns=self.velocities.columns, header=False))        
        cfile.write("\n")


        #Bonds
        cfile.write('\nBonds\n\n')
        cfile.write(self.bonds.to_string(index=False, columns=self.bonds.columns, header=False))
        cfile.write("\n")

        
        #Angles
        cfile.write('\nAngles\n\n') 
        cfile.write(self.angles.to_string(index=False, columns=self.angles.columns, header=False))
        cfile.write("\n")


        #Dihedrals
        cfile.write('\nDihedrals\n\n') 
        cfile.write(self.dihedrals.to_string(index=False, columns=self.dihedrals.columns, header=False))
        cfile.write("\n")


        #Impropers
        cfile.write('\nImpropers\n\n') 
        cfile.write(self.impropers.to_string(index=False, columns=self.impropers.columns, header=False))
        cfile.write("\n")



if __name__ == "__main__":  # tests
    from NAMDdata import NAMDdata
    l = LammpsData()
    ch = NAMDdata()
    ch.readFiles("2rvd_autopsf.pdb", "2rvd_autopsf.psf", "par_all36_prot.prm")

    l.loadNAMDdata(ch)
    l.writeConf("2rvd.data")
