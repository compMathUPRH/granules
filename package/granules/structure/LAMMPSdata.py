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

#try:
#    import networkx as nx
#except Exception as e:
#    print("module networkx not found, some functionality may not be available")

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

#Clase para generar el archivo de configuracion para simulaciones en Lammps
class InFileGenerator():
    '''Write the default configuration for a .in file'''
    Commentarios = ['Generar Caja con Atomos','Especificaciones de la corrida','Vecinos','Output y Simulacion']
    
    LammpsConfDefault = {'variable':'t index 5000','units':'lj','atom_style':'atomic','lattice':'fcc 1.0',
                         'region':'box block 0  10.0 0 10.0 0 10.0','create_box':'3 box',
                         'create_atoms':['1 single 5 5 8','1 single 8 5 6','2 single 7 5 7','2 single 3 6 6',
                         '3 single 3 5 3','3 single 9 5 7'],'mass':['1 460','2 70','3 10'],
                         'velocity':'all create 4 87287 loop geom','pair_style':'lj/cut 2.5',
                         'pair_coeff':['1 1 5.0 2.0 10.0','1 2 7.5 1.5 10.0','2 2 10.0 1.0 10.0',
                         '2 3 3.0 1.2 10.0','3 3 6.0 2.0 10.0'],'neighbor':'0.3 bin',
                         'neigh_modify':'delay 0 every 40 check no','fix':'1 all nve','thermo':'200',
                         'thermo_style':'custom elapsed pe ke etotal press temp',
                         'dump':'d0 all image 100 dump.\*.jpg type type','run':'$t'}
    
    def __init__(self,filename='in.Prueba1'):
        self.filename = filename
            
    def writeFile(self):
        '''Write configurations in txt file'''
        file = open(self.filename,"w+")
        cont = 0
        for variables in self.LammpsConfDefault:
            if type(self.LammpsConfDefault[variables]) != list:
                if variables not in ['lattice','velocity','neighbor','fix']: 
                    file.write(variables +' '+self.LammpsConfDefault[variables]+'\n')
                    print(variables +' '+self.LammpsConfDefault[variables]+'\n')
                else:
                    file.write('\n########'+self.Commentarios[cont]+'#######\n')
                    file.write(variables +' '+self.LammpsConfDefault[variables]+'\n')
                    print('\n'+variables +' '+self.LammpsConfDefault[variables]+'\n')
                    cont+=1
                    
            else:
                print("#########Key con listas###########")
                for elemt in self.LammpsConfDefault[variables]:
                    file.write(variables +' '+elemt+'\n')
                    print(variables +' '+elemt+'\n')
                print("\n")
        file.close()
        print("Termine de escribir.")
        #return file
        
#===================================================================
#Clases temporeras que hereda de panda,las subclases de sus clases originales(nombre)
#   se pasaron a estas nuevas clases.
"""
class Atoms():#cambiado LammpsDataFrame
    def __init__(self):
        pass
   
class MolecularTopology():
    def __init__(self):
        pass
    
class ForceField():
    def __init__(self):
         pass
"""
#===================================================================
         
class LammpsDataFrame(pd.DataFrame):
    '''
    @property
    def _constructor(self):
        return LammpsDataFrame
    '''
    
    @property
    def _constructor(self):
        '''In absence of this method many operations would return a DataFrame, not an type(self)
        '''
        #print("AtomsDF._constructor", type(self))
        return type(self)
        #return AtomsFull.atomsDF


    def __getitem__(self, key):
        ''' Selecting columns no longer produces an AtomDF. It will return a
        DataFrame or Slice'''
        return pd.DataFrame(data=self)[key]


    def add(self, data):
        print("CUIDADO: Llegó a LammpsDataFrame.add(). Esta tiene que ser sustituida.")
        raise Exception("CUIDADO: Llegó a LammpsDataFrame.add(). Esta tiene que ser sustituida.")
        


#===================================================================

class AtomsFull:
    def __init__(self):
        
        # atom-property sections
        self.atoms       = self.AtomsDF()
        self.velocities  = self.VelocitiesDF()
        self.masses      = self.MassesDF()
               
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

    class AtomsDF(LammpsDataFrame):
        ''' This subclass of pd.DataFrame will always have the columns specified
            at method __init__. Slicing this table will return DataFrame or Slice.
        '''
        def __init__(self, *args, **kwargs):
            column_names = ['aID', 'Mol_ID', 'aType', 'Q', 'x', 'y', 'z', 'Nx', 'Ny', 'Nz']
            super(AtomsFull.AtomsDF,self).__init__(*args, **dict(kwargs, columns=column_names))

            #print("AtomsDF len(args) ", len(args), kwargs.keys())
            #if len(args) > 0: print("AtomsDF type(args[0])=", type(args[0]))
            if not hasattr(self, "atomTypes"): 
                self.atomTypes = {}
                
            if not hasattr(self, "moleculeTypes"): 
                self.moleculeTypes = {}

            if len(args) > 0 and type(args[0]) == type(self):
                self.atomTypes = args[0].atomTypes
                self.moleculeTypes = args[0].moleculeTypes

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
            sel_pdb     = charmm.pdb[['ID','x','y','z']].set_index('ID')
            sel         = sel_pdb.join(sel_psf)
            sel['aType'] = sel_psf['Type'].map(charmmTypeToInt)
            sel.reset_index(inplace=True)
            sel         .rename(columns={"Charge":"Q",'ID':'aID'}, inplace=True)
    
            # add remining columns
            if hasattr(charmm, 'atomsInMolecules'):
                self.moleculeTypes = {v:k+1 for k,v in enumerate(set(charmm.atomsInMolecules.molname))}
                sel = sel.join(charmm.atomsInMolecules.set_index('ID'), on='aID')
                sel['Mol_ID'] = sel['molname'].map(self.moleculeTypes)
            else:
                sel['Mol_ID'] = np.ones((len(sel), 1), dtype=np.int8)
            sel['Nx']     = np.zeros((len(sel), 1))
            sel['Ny']     = np.zeros((len(sel), 1))
            sel['Nz']     = np.zeros((len(sel), 1))
            #sel['ID']     = np.arange(1, len(sel)+1)

            # rearrange columns
            sel = sel[['aID', 'Mol_ID', 'aType', 'Q', 'x', 'y', 'z', 'Nx', 'Ny', 'Nz']]
            
            #sel.reset_index(inplace=True)
            self.__init__(sel)

            self.astype({
                         'Mol_ID' :int,
                         'aType' :int,
                         'Q' :float,
                         'x' :float,
                         'y' :float,
                         'z' :float,
                         'Nx' :int,
                         'Ny' :int,
                         'Nz' : int
                        })
            self.atomTypes = {v: k for k, v in charmmTypeToInt.items()}
       
        def center(self):
            return self[['x', 'y', 'z']].mean()
    
        def move(self, b):
            for coord in ['x', 'y', 'z']:
                self[coord] = self[coord] + b[coord]
            return self
    
        def update_from_file(self,archivo):
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
        
        def update_corrdinates(self, coords):
            ''' Update coordintes. Coordinates in self may be less that the ones 
            received in coords because a selection may have been made.
            
            coords: tuple with Series for X, Y, and Z.
            '''
            #print(self)
            data = pd.DataFrame(np.array(coords).T,columns = ['x', 'y', 'z'])
            data.set_index(np.arange(1,data.shape[0]+1), inplace=True)
            data = data[data.index.isin(self.aID)]
            for column in ['x', 'y', 'z']: self[column] = data[column].values
            #print(data)
         
    class MassesDF(LammpsDataFrame):
        ''' This subclass of pd.DataFrame will always have the columns specified
            at method __init__. Slicing this table will return DataFrame or Slice.
        '''
        def __init__(self, *args, **kwargs):
            column_names = ['aType', 'Mass']
            super(AtomsFull.MassesDF,self).__init__(*args, **dict(kwargs, columns=column_names))
        
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
            sel      = sel_self.join(sel_psf).drop_duplicates()
    
            # rename columns
            #sel      .drop(columns='aID', inplace=True)
            sel.reset_index(drop=True)
            #print(sel.dtypes)
    
            self.__init__(sel)

    
   
    class VelocitiesDF(LammpsDataFrame):
        ''' This subclass of pd.DataFrame will always have the columns specified
            at method __init__. Slicing this table will return DataFrame or Slice.
        '''
        def __init__(self, *args, **kwargs):
            column_names = ['vID', 'Vx', 'Vy', 'Vz']
            super(AtomsFull.VelocitiesDF,self).__init__(*args, **dict(kwargs, columns=column_names))


        def setToZero(self, atoms):
            ''' Sets velocities of atoms to zero.
    
            Parameter
            -----------------
            atoms     : AtomsDF
                an AtomsDF object
            '''
    
            # extract info from LAMMPS
            #print("VelocitiesDF.setToZero", type(atoms), atoms)
            sel = atoms[['aID']].copy().rename(columns={'aID':'vID'})
            #sel.rename(columns={'aID':'vID'}, inplace=True)
            sel['Vx']     = np.zeros((len(sel), 1))
            sel['Vy']     = np.zeros((len(sel), 1))
            sel['Vz']     = np.zeros((len(sel), 1))
            #print("VelocitiesDF sel = ", sel.dtypes)
    
            self.__init__(sel)

 
 
class MolecularTopologyData:
    def __init__(self):
        
        # molecular topology sections
        self.angles      = self.AnglesDF()
        self.bonds       = self.BondsDF()
        self.dihedrals   = self.DihedralsDF()
        self.impropers   = self.ImpropersDF()
        
    def setFromNAMD(self,charmm): 
        '''Llama a la funcion setFromNAMD() de las clases de la clase MolecularTopolyData,
            asignadas en los atributo.'''
        
        # molecular topology sections
        self.bonds.setFromNAMD(charmm)
        self.angles.setFromNAMD(charmm)
        self.dihedrals.setFromNAMD(charmm)
        self.impropers.setFromNAMD(charmm)

    class AnglesDF(LammpsDataFrame):
        ''' This subclass of pd.DataFrame will always have the columns specified
            at method __init__. Slicing this table will return DataFrame or Slice.
        '''
        def __init__(self, *args, **kwargs):
            column_names = ['anID', 'anType', 'atom1', 'atom2', 'atom3']
            super(MolecularTopologyData.AnglesDF,self).__init__(*args, **dict(kwargs, columns=column_names))


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
            angles['atom1'] = charmm.psf.angles.copy()['atom1']
            angles['atom2'] = charmm.psf.angles.copy()['atom2']
            angles['atom3'] = charmm.psf.angles.copy()['atom3']
    
            angles.rename(columns={'ID':'anID', 'Type':'anType'}, inplace=True)
            #print(angles)
    
            super(MolecularTopologyData.AnglesDF, self).__init__(angles)
            #print(self.dtypes)
    
    class BondsDF(LammpsDataFrame):
        ''' This subclass of pd.DataFrame will always have the columns specified
            at method __init__. Slicing this table will return DataFrame or Slice.
        '''
        def __init__(self, *args, **kwargs):
            column_names = ['bID', 'bType', 'atom1', 'atom2']
            super(MolecularTopologyData.BondsDF,self).__init__(*args, **dict(kwargs, columns=column_names))


        def selectTypes(self, bondType):
            return self[self.bType == bondType].bID
    
    
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
            bonds['atom1'] = charmm.psf.bonds.copy()['atom1']
            bonds['atom2'] = charmm.psf.bonds.copy()['atom2']
    
            bonds.rename(columns={'ID':'bID', 'Type':'bType'}, inplace=True)
    
            super(MolecularTopologyData.BondsDF, self).__init__(bonds)
            
    
           
    class DihedralsDF(LammpsDataFrame):
        ''' This subclass of pd.DataFrame will always have the columns specified
            at method __init__. Slicing this table will return DataFrame or Slice.
        '''
        def __init__(self, *args, **kwargs):
            column_names = ['dID', 'dType', 'atom1', 'atom2', 'atom3', 'atom4']
            super(MolecularTopologyData.DihedralsDF,self).__init__(*args, **dict(kwargs, columns=column_names))

        
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
            dihes['atom1'] = charmm.psf.dihedrals.copy()['atom1']
            dihes['atom2'] = charmm.psf.dihedrals.copy()['atom2']
            dihes['atom3'] = charmm.psf.dihedrals.copy()['atom3']
            dihes['atom4'] = charmm.psf.dihedrals.copy()['atom4']
    
            dihes.rename(columns={'ID':'dID', 'Type':'dType'}, inplace=True)
            #print(dihes)
    
            super(MolecularTopologyData.DihedralsDF, self).__init__(dihes)
            #print(self.dtypes)
    
    
    class ImpropersDF(LammpsDataFrame):
        ''' This subclass of pd.DataFrame will always have the columns specified
            at method __init__. Slicing this table will return DataFrame or Slice.
        '''
        def __init__(self, *args, **kwargs):
            column_names = ['iID', 'iType', 'atom1', 'atom2', 'atom3', 'atom4']
            super(MolecularTopologyData.ImpropersDF,self).__init__(*args, **dict(kwargs, columns=column_names))

        
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
            impros['atom1'] = charmm.psf.impropers.copy()['atom1']
            impros['atom2'] = charmm.psf.impropers.copy()['atom2']
            impros['atom3'] = charmm.psf.impropers.copy()['atom3']
            impros['atom4'] = charmm.psf.impropers.copy()['atom4']
    
            impros.rename(columns={'ID':'iID', 'Type':'iType'}, inplace=True)
            #print(impros)
    
            super(MolecularTopologyData.ImpropersDF, self).__init__(impros)
            #print(self.dtypes)



class CharmmForceField:
     def __init__(self):

        #force field sections
        #self.angleAngleCoeffs       
        #self.angleAngleTorsionCoeffs
        #self.angleTorsionCoeffs     
	
        #self.bondAngleCoeffs        
        #self.bondBond13Coeffs       
        
        #self.middleBondTorsion      
        #self.endBondTorsionCoeffs  
        
        #Moverlos a LAmmpsData()
        self.angleCoeffs           = self.AngleCoeffs()
        self.bondCoeffs             = self.BondCoeffs()
        
        self.dihedralCoeffs         = self.DihedralCoeffs()
        self.improperCoeffs         = self.ImproperCoeffs()
        self.pairCoeffs             = self.PairCoeffs()
        
        
     def setFromNAMD(self,charmm,atompropertydata,topology): #añadi este codigo nuevo 
        
        self.pairCoeffs.setFromNAMD(charmm, atompropertydata.atoms)
        self.bondCoeffs.setFromNAMD(charmm, topology.bonds)
        self.angleCoeffs.setFromNAMD(charmm, topology.angles)
        self.dihedralCoeffs.setFromNAMD(charmm, topology.dihedrals)
        self.improperCoeffs.setFromNAMD(charmm, topology.impropers)
        
     def NonBondEnergy(self,atompropertydata,topology):
        ''' Computes CHARMM Lennard-Jones energy.
            Formula: Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
                    Eps,i,j = sqrt(eps,i * eps,j)
                    Rmin,i,j = Rmin/2,i + Rmin/2,j

            Computes CHARMM Coulumb energy.
            Formula: qi qj/rij / epsilon_0

            returns (L-J, Coulomb)
        '''
        #from scipy.constants import epsilon_0, physical_constants

        NONB_CUTOFF = 13.0

        print("CharmmForceField.NonBondEnergy")
        
        # generate all pairs of atoms IDs
        atoms = atompropertydata.atoms.copy() #Cambio
        atomIds = atoms[['aID']]
        print("CharmmForceField.NonBondEnergy atomIds.columns=", atomIds.columns)
        atomIds = atomIds.assign(key=np.ones(len(atomIds)))
        #atomIds['key'] = np.ones(len(atomIds))
        #print("CharmmForceField.NonBondEnergy atomIds.columns=", atomIds.columns)
        atomIds = pd.merge(atomIds, atomIds, on='key')[['aID_x', 'aID_y']]
        atomIds = atomIds[atomIds['aID_x'] < atomIds['aID_y']]

        atomIds['nbID'] = np.arange(len(atomIds))

        # compute pairwise distances
        print("CharmmForceField.NonBondEnergy len(atomIds) =", len(atomIds))
        from scipy.spatial.distance import pdist
        atomIds['rij'] = pdist(atoms.set_index('aID')[['x', 'y', 'z']].values)
        atomIds = atomIds[atomIds['rij'] < NONB_CUTOFF]

        # remove bonded atoms
        atomIds['p'] = list(zip(atomIds.aID_x, atomIds.aID_y))
        bonds = topology.bonds.copy()
        bonds['p'] = list(zip(bonds.atom1, bonds.atom2))
        atomIds = atomIds.set_index('p').join(bonds.set_index('p'))
        atomIds = atomIds[atomIds.bID.isna()][['aID_x','aID_y','nbID','rij']].reset_index(drop=True)
        del bonds

        # remove angled atoms
        atomIds['p'] = list(zip(atomIds.aID_x, atomIds.aID_y))
        angles = topology.angles.copy()
        angles['p'] = list(zip(angles.atom1, angles.atom3))
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
        print("CharmmForceField.NonBondEnergy len(atomIds (clean)) =", len(atomIds))

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
        print("CharmmForceField.NonBondEnergy END")
        #print(atomIds)
        # return LJ and Coulomb
        COULOMB = 332.0636

        '''
        return -np.sum(atomIds.epsilon * (atomIds.rij**12 - atomIds.rij**6)), \
               np.sum((physical_constants['electric constant'][0] * atomIds.qi * atomIds.qj) / (atomIds.rij * epsilon_0))
        '''
        return -np.sum(atomIds.epsilon * ((atomIds.sigma/atomIds.rij)**12 - (atomIds.sigma/atomIds.rij)**6)), \
               COULOMB * np.sum((atomIds.qi * atomIds.qj) / (atomIds.rij))

     def BondEnergy(self,atompropertydata,topology):
        ''' Computes CHARMM bond energy.

            Formula: sum K * (bij - b0)**2
        '''
        bi = topology.bonds.set_index('atom1').join(
                atompropertydata.atoms.set_index('aID')
             )[['bID','x', 'y', 'z']].set_index('bID')

        bj = topology.bonds.set_index('atom2').join(
                atompropertydata.atoms.set_index('aID')
             )[['bID','x', 'y', 'z']].set_index('bID')


        bij = np.sqrt(np.sum((bi - bj) ** 2, axis='columns'))  # bonds lengths

        coeffs = topology.bonds[['bID','bType']].set_index('bType').join(
                        self.bondCoeffs.set_index('bType')
                 ).reset_index(drop=True)
        K = coeffs[['bID','Spring_Constant']].set_index('bID').Spring_Constant
        b0 = coeffs[['bID','Eq_Length']].set_index('bID').Eq_Length

        return np.sum(K * (bij-b0)**2)


     def AngleEnergy(self,atompropertydata,topology):
        ''' Computes CHARMM angle energy.
            Formula: sum K * (aij - a0)**2

        '''
        bi = topology.angles.set_index('atom1').join(
                atompropertydata.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        bj = topology.angles.set_index('atom2').join(
                atompropertydata.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        bk = topology.angles.set_index('atom3').join(
                atompropertydata.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        # compute angles
        l1 = bi - bj
        l2 = bk - bj

        norm1 = np.sqrt(np.square(l1).sum(axis=1))
        norm2 = np.sqrt(np.square(l1).sum(axis=1))
        dot = l1.x * l2.x + l1.y * l2.y + l1.z * l2.z
        angles = np.arccos(dot / (norm1 * norm2))

        coeffs = topology.angles[['anID','anType']].set_index('anType').join(
                        self.angleCoeffs.set_index('anType')
                 ).reset_index(drop=True)

        K = coeffs[['anID','Ktheta']].set_index('anID').Ktheta
        a0 = np.radians(coeffs[['anID','Theta0']].set_index('anID').Theta0)

        return np.sum(K * (angles-a0)**2)

     def DihedralsEnergy(self,atompropertydata,topology):
        ''' Computes CHARMM angle energy.
            Formula: sum K * (1 + cos(n * x - d))
        '''
        bi = topology.dihedrals.set_index('atom1').join(
                atompropertydata.atoms.set_index('aID')
             )[['dID','x', 'y', 'z']].set_index('dID')

        bj = topology.dihedrals.set_index('atom2').join(
                atompropertydata.atoms.set_index('aID')
             )[['dID','x', 'y', 'z']].set_index('dID')

        bk = topology.dihedrals.set_index('atom3').join(
                atompropertydata.atoms.set_index('aID')
             )[['dID','x', 'y', 'z']].set_index('dID')

        bl = topology.dihedrals.set_index('atom4').join(
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
        print(topology.dihedrals)
        coeffs = topology.dihedrals[['dID','dType']].set_index('dType').join(
                        self.dihedralCoeffs.set_index('dType')
                 ).reset_index(drop=True)
        K = coeffs[['dID','Kchi']].set_index('dID').Kchi
        n = coeffs[['dID','n']].set_index('dID').n
        d = coeffs[['dID','delta']].set_index('dID').delta



        return np.sum(K * (1 + np.cos(n * angles - d)))
        '''
        '''


     def NonBondForce(self,atompropertydata,topology):
        ''' Computes CHARMM Lennard-Jones energy.
            Formula: Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
                    Eps,i,j = sqrt(eps,i * eps,j)
                    Rmin,i,j = Rmin/2,i + Rmin/2,j

            Computes CHARMM Coulumb energy.
            Formula: qi qj/rij / epsilon_0

            returns (L-J, Coulomb)
        '''
        NONB_CUTOFF = 13.0
        print("CharmmForceField.NonBondForce()")
        
        # generate all pairs of atoms IDs
        #atomIds = atompropertydata.atoms[['aID']].copy()#Fastidia a pyreverse
        raise Exception("Instruction in CharmmForceField.NonBondForce() hidden to make pyreverse work.")
        
        atomIds['key'] = np.ones(len(atomIds))
        atomIds = pd.merge(atomIds, atomIds, on='key')[['aID_x', 'aID_y']]
        atomIds = atomIds[atomIds['aID_x'] < atomIds['aID_y']]

        atomIds['nbID'] = np.arange(len(atomIds))

        # compute pairwise distances
        print('CharmmForceField.NonBondForce: len(atomIds)=', len(atomIds))
        from scipy.spatial.distance import pdist
        atomIds['rij'] = pdist(atompropertydata.atoms.set_index('aID')[['x', 'y', 'z']].values)
        atomIds = atomIds[atomIds['rij'] < NONB_CUTOFF]

        # remove bonded atoms
        print('CharmmForceField.NonBondForce: len(atomIds < NONB_CUTOFF)=', len(atomIds))
        atomIds['p'] = list(zip(atomIds.aID_x, atomIds.aID_y))
        bonds = topology.bonds.copy()
        bonds['p'] = list(zip(bonds.atom1, bonds.atom2))
        atomIds = atomIds.set_index('p').join(bonds.set_index('p'))
        atomIds = atomIds[atomIds.bID.isna()][['aID_x','aID_y','nbID','rij']].reset_index(drop=True)
        del bonds
        print('CharmmForceField.NonBondForce: len(atomIds < NONB_CUTOFF -  BONDS)=', len(atomIds))

        # remove angled atoms
        atomIds['p'] = list(zip(atomIds.aID_x, atomIds.aID_y))
        angles = topology.angles.copy()
        angles['p'] = list(zip(angles.atom1, angles.atom3))
        atomIds = atomIds.set_index('p').join(angles.set_index('p'))
        atomIds = atomIds[atomIds.anID.isna()][['aID_x','aID_y','nbID','rij']].reset_index(drop=True)
        del angles
        print('CharmmForceField.NonBondForce: len(atomIds < NONB_CUTOFF -  BONDS - ABGLES)=',len(atomIds))

        # get atom types and charges
        atomIds = atomIds.set_index('aID_x', drop=False).join(atompropertydata.atoms[['aID', 'Q']].set_index('aID'))
        atomIds.rename(columns={'Q':'qi'}, inplace=True)
        atomIds = atomIds.join(atompropertydata.atoms[['aID', 'aType']].set_index('aID')).reset_index(drop=True)
        atomIds.rename(columns={'aType':'aiType'}, inplace=True)
        atomIds = atomIds.set_index('aID_y', drop=False).join(atompropertydata.atoms[['aID', 'Q']].set_index('aID'))
        atomIds = atomIds.join(atompropertydata.atoms[['aID', 'aType']].set_index('aID')).reset_index(drop=True)
        atomIds.rename(columns={'aType':'ajType', 'Q':'qj'}, inplace=True)

        # get epsilons and sigmas for each atom type
        #areglar error
        sameTypes = self.pairCoeffs[ self.pairCoeffs['aType'] == self.pairCoeffs['aType2'] ]
        sameTypes.drop('aType2', axis=1, inplace=True) 
        
        atomIds = atomIds.set_index('aiType').join(
                        sameTypes.set_index('aType')
                 ).reset_index(drop=True)
        atomIds.drop(columns=['epsilon1_4', 'sigma1_4'], inplace=True)
        atomIds.rename(columns={'epsilon':'epsilon_i', 'sigma':'sigma_i'}, inplace=True)
        atomIds = atomIds.set_index('ajType').join(
                        sameTypes.set_index('aType')
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


     def BondForce(self,atompropertydata,topology):
        ''' Computes CHARMM bond energy.
            Formula: sum K * (bij - b0)**2

        '''
        bi = topology.bonds.set_index('atom1').join(
                atompropertydata.atoms.set_index('aID')
             )[['bID','x', 'y', 'z']].set_index('bID')

        bj = topology.bonds.set_index('atom2').join(
                atompropertydata.atoms.set_index('aID')
             )[['bID','x', 'y', 'z']].set_index('bID')

        bij = bj - bi  # bonds
        rbij = np.sqrt(np.sum((bij) ** 2, axis='columns'))  # bonds lengths

        coeffs = topology.bonds[['bID','bType']].set_index('bType').join(
                        self.bondCoeffs.set_index('bType')
                 ).reset_index(drop=True)
        K = coeffs[['bID','Spring_Constant']].set_index('bID').Spring_Constant
        b0 = coeffs[['bID','Eq_Length']].set_index('bID').Eq_Length

        forces = -2 * K * (rbij-b0)
        
        bij = bij.div(rbij, axis=0)  # normalize
        wi = bij.mul(forces, axis=0)
        #wj = bij.mul(-forces, axis=0)
        
        fi = wi.join(topology.bonds.set_index('bID'))[['x','y','z','atom1']].groupby('atom1').sum()
        fj = wi.join(topology.bonds.set_index('bID'))[['x','y','z','atom2']].groupby('atom2').sum()
        fi.index.names = ['aID']
        fj.index.names = ['aID']
        ff = fi.add(fj, axis=0, fill_value=0)

        return ff



     def AngleForce(self,atompropertydata,topology):
        ''' Computes CHARMM angle energy.
            Formula: sum K * (aij - a0)**2
        '''
        
        print("CharmmForceField.AngleForce")
        bi = topology.angles.set_index('atom1').join(
                atompropertydata.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        bj = topology.angles.set_index('atom2').join(
                atompropertydata.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        bk = topology.angles.set_index('atom3').join(
                atompropertydata.atoms.set_index('aID')
             )[['anID','x', 'y', 'z']].set_index('anID')

        # compute angles
        l1 = bi - bj
        l2 = bk - bj

        norm1 = np.sqrt(np.square(l1).sum(axis=1))
        norm2 = np.sqrt(np.square(l1).sum(axis=1))
        dot = l1.x * l2.x + l1.y * l2.y + l1.z * l2.z
        angles = np.arccos(dot / (norm1 * norm2))

        coeffs = topology.angles[['anID','anType']].set_index('anType').join(
                        self.angleCoeffs.set_index('anType')
                 ).reset_index(drop=True)

        K = coeffs[['anID','Ktheta']].set_index('anID').Ktheta
        a0 = np.radians(coeffs[['anID','Theta0']].set_index('anID').Theta0)

        c11 = l1.x * l1.x + l1.y * l1.y + l1.z * l1.z
        c12 = l1.x * l2.x + l1.y * l2.y + l1.z * l2.z
        c22 = l2.x * l2.x + l2.y * l2.y + l2.z * l2.z
        cd = np.sqrt(c11 * c22)
        #c = c12 / cd
        
        forces = -2 * K * (angles-a0)
        
        w1 = l1.mul(c12 / c11, axis=0).add(-l2, axis=0).mul(forces / cd, axis=0)
        w2 = l1.sub(l2.mul(c12 / c22, axis=0),  axis=0).mul(forces / cd, axis=0)

        #print(self.angles.columns.values)
        f1 = w1.join(topology.angles.set_index('anID'))[['x','y','z','atom1']].groupby('atom1').sum()
        f2 = w1.sub(w2, axis=0).join(topology.angles.set_index('anID'))[['x','y','z','atom2']].groupby('atom2').sum()
        f3 = w2.join(topology.angles.set_index('anID'))[['x','y','z','atom3']].groupby('atom3').sum()
        f1.index.names = ['aID']
        f2.index.names = ['aID']
        f3.index.names = ['aID']

        ff = f1.add(f2.add(f3, axis=0, fill_value=0), axis=0, fill_value=0)
        #print (ff)
        print("CharmmForceField.AngleForce  END")
        return ff

     def Force(self,atompropertydata,topology):
        print("CharmmForceField.Force()")
        return self.NonBondForce(atompropertydata,topology).add(
                self.BondForce(atompropertydata,topology), axis=0).add(
                self.AngleForce(atompropertydata,topology), axis=0)
        
     def Energy(self,atompropertydata,topology):
        return self.NonBondEnergy(atompropertydata,topology).add(
                self.BondEnergy(atompropertydata,topology)).add(
                self.AngleEnergy(atompropertydata,topology)).add(
                self.DihedralsEnergy(atompropertydata,topology))

     class PairCoeffs(LammpsDataFrame):
         def __init__(self,data=None, dtype=None, copy=False):
            super(CharmmForceField.PairCoeffs, self).__init__(data=data, columns=['aType1','aType2', 'epsilon', 'sigma', 'epsilon1_4', 'sigma1_4'], dtype=dtype, copy=copy)
    
         def setFromNAMD(self, charmm, atoms):
            ''' Extracts info from PRM and PSF objects into self.
    
            Parameter
            -----------------
            charmm : NAMDdata
                NAMDdata object
    
            mass : MassDF
                AtomsDF object associateed with these PairCoeffs
            '''
            from itertools import product

            # extract info from charmm
            psf_types   = charmm.psf.atoms[['ID', 'Type']].set_index('ID').to_dict()['Type']
            #print("PairCoeffs psf_types ", psf_types)
    
            # substitute atoms numbers with charmm atom types
            nonbonded       = atoms[['aID', 'aType']].copy()

            nonbonded['types'] = nonbonded.aID.map(psf_types)
            #print("PairCoeffs nonbonded B \n", nonbonded)
            nonbonded.drop(columns=['aID'], inplace=True)
            nonbonded.drop_duplicates(inplace=True)
            #print("PairCoeffs nonbonded C \n", nonbonded)
    
            prmFF = charmm.prm.nonbonded.getCoeffs()
            #print("Aqui el DF: ",prmFF)
    
            # add charge and energy to atoms
            nonbonded['epsilon'] = nonbonded.types.map(prmFF.epsilon.to_dict())
            nonbonded['sigma'] = nonbonded.types.map(prmFF.Rmin2.to_dict())
            nonbonded['epsilon1_4'] = nonbonded.types.map(prmFF.epsilon.to_dict())
            nonbonded['sigma1_4'] = nonbonded.types.map(prmFF.Rmin2.to_dict())
            nonbonded.drop(columns=['types'], inplace=True)
            
            nonbonded = nonbonded[['aType', 'epsilon','sigma', 'epsilon1_4', 'sigma1_4']].set_index('aType', drop=True)
            #print("PairCoeffs aqui el nonbonded, la tabla:\n",nonbonded)
            
            # generate all pairs
            pairs = pd.Series([i for i in product(nonbonded.index, repeat=2) if i[0] <= i[1]])
            pairs = pd.concat([pairs.apply(lambda x: x[0]), pairs.apply(lambda x: x[1])],1).rename(columns={0:'aType1', 1:'aType2'})
            pairs['epsilon'] = -np.sqrt(pairs.aType1.map(nonbonded.epsilon) * pairs.aType2.map(nonbonded.epsilon))
            pairs['epsilon1_4'] = -np.sqrt(pairs.aType1.map(nonbonded.epsilon1_4) * pairs.aType2.map(nonbonded.epsilon1_4))
            pairs['sigma'] = (pairs.aType1.map(nonbonded.sigma) + pairs.aType2.map(nonbonded.sigma))/2
            pairs['sigma1_4'] = (pairs.aType1.map(nonbonded.sigma1_4) + pairs.aType2.map(nonbonded.sigma1_4))/2
            #print("\nPairCoeffs pairs:\n",pairs)
            self.__init__(pairs)
            #print("\nPairCoeffs Nans:\n",nonbonded.isna().sum())
            #print(self.dtypes

         def dropPair(self, atomType1, atomType2):
             ''' Removes the potential parameters for a pair (atomType1, atomType2).
             '''
             index_names = self[(self['aType1'] == atomType1) &   (self['aType2'] == atomType2)].index
             self.drop(index_names, inplace = True) 

         def setPair(self, atomType1, atomType2, epsilon, sigma, epsilon1_4, sigma1_4):
             ''' sets the potential parameters for a pair (atomType1, atomType2).
                 If the pair exist it will replace it.
             '''
             if atomType1 > atomType2: atomType1, atomType2 = atomType2, atomType1
             
             self.dropPair(atomType1, atomType2) 
             self.__init__(
                     self.append({'aType1':atomType1,'aType2':atomType2,'epsilon':epsilon,'sigma':sigma,'epsilon1_4':epsilon1_4,'sigma1_4':sigma1_4}, sort=False, ignore_index=True)
                  )   

     class AngleCoeffs(LammpsDataFrame):
         def __init__(self,data=None, dtype=None, copy=False):
             super(CharmmForceField.AngleCoeffs, self).__init__(data =data, columns=['anType', 'Ktheta', 'Theta0', 'Kub', 'S0'], dtype=dtype, copy=copy)
    
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
            angles     = angles[['anType', 'atom1', 'atom2', 'atom3']]
            angles.atom1 = angles.atom1.map(psf_types)
            angles.atom2 = angles.atom2.map(psf_types)
            angles.atom3 = angles.atom3.map(psf_types)
            #print("CharmmForceField.AngleCoeffs{}\n{}".format(type(angles), angles))
            angles['atuple'] = list(zip(angles.atom1, angles.atom2, angles.atom3))
            #print("CharmmForceField.AngleCoeffs{}\n{}".format(type(angles), angles))
            angles.drop(columns=['atom1', 'atom2', 'atom3'], inplace=True)
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
    
            self.__init__(angles)
            #print("\nAngleCoeffs Nans:\n",angles.isna().sum())
            #print("Nans:\n",self)
    
    
    
     class BondCoeffs(LammpsDataFrame):
         def __init__(self,data=None, dtype=None, copy=False):
             super(CharmmForceField.BondCoeffs, self).__init__(data=data, columns=['bType','Spring_Constant','Eq_Length'], dtype=dtype, copy=copy)
    
         def setFromNAMD(self, charmm, bondsPar):
            ''' Extracts info from PRM and PSF objects into self.
    
            Parameter
            -----------------
            charmm : NAMDdata
                NAMDdata object
    
            bonds : BondsDF
                BondsDF object associateed with these BondCoeffs
            '''
    
            # extract info from charmm
            #print("\nBondCoeffs bondsPar:\n",bondsPar)
            psf_types   = charmm.psf.atoms[['ID', 'Type']].set_index('ID').to_dict()['Type']
            #print("\nBondCoeffs psf_types:\n",psf_types)
    
            # substitute atoms numbers with charmm atom types
            #bonds     = bondsPar.copy()
            bonds     = bondsPar[['bID', 'bType', 'atom1', 'atom2']]
            
            bonds.atom1 = bonds.atom1.map(psf_types)
            bonds.atom2 = bonds.atom2.map(psf_types)
            bonds['atuple'] = list(zip(bonds.atom1, bonds.atom2))
            #print("\nBondCoeffs bonds:\n",bonds)
            bonds.drop(columns=['bID', 'atom1', 'atom2'], inplace=True)
            bonds.drop_duplicates(inplace=True)
            #print(bonds)
    
            prmFF = charmm.prm.bonds.getCoeffs()
            #print("\nBondCoeffs prmFF:\n",prmFF)
    
    
            # add Kb and b0 to bonds
            bonds['Spring_Constant'] = bonds.atuple.map(prmFF.Kb.to_dict())
            bonds['Eq_Length'] = bonds.atuple.map(prmFF.b0.to_dict())
            bonds.drop(columns=['atuple'], inplace=True)
            #print("\nBondCoeffs bonds:\n",bonds)
            #print("\nBondCoeffs Nans:\n",bonds.isna().sum())
    
            bonds.rename(columns={'bID':'bType'}, inplace=True)
    
            self.__init__(bonds)
            #print(self.dtypes)
    
    
     class DihedralCoeffs(LammpsDataFrame):
         def __init__(self,data=None, dtype=None, copy=False):
             super(CharmmForceField.DihedralCoeffs, self).__init__(data=data, columns=['dType', 'Kchi', 'n', 'delta'], dtype=dtype, copy=copy)
    
         def setFromNAMD(self, charmm, dihedrals_par):
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
    
            #dihedrals     = dihedrals_par.copy()
            dihedrals     = dihedrals_par[['atom1', 'atom2', 'atom3', 'atom4']]
            
            dihedrals.atom1 = dihedrals.atom1.map(psf_types)
            dihedrals.atom2 = dihedrals.atom2.map(psf_types)
            dihedrals.atom3 = dihedrals.atom3.map(psf_types)
            dihedrals.atom4 = dihedrals.atom4.map(psf_types)
            dihedrals['atuple'] = list(zip(dihedrals.atom1, dihedrals.atom2, dihedrals.atom3, dihedrals.atom4))
            dihedrals.drop(columns=['atom1', 'atom2', 'atom3', 'atom4'], inplace=True)
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
            
            self.__init__(dihedrals.astype({'Kchi':float, 'n':int, 'delta':int,'Weighting_Factor':float}))
    
    
     class ImproperCoeffs(LammpsDataFrame):
         def __init__(self,data=None, dtype=None, copy=False):
             super(CharmmForceField.ImproperCoeffs, self).__init__(data=data, columns=['iType', 'Kchi', 'n', 'delta'], dtype=dtype, copy=copy)
    
         def setFromNAMD(self, charmm, impropers_par):
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
    
            #impropers     = impropers_par.copy()
            impropers     = impropers_par[['atom1', 'atom2', 'atom3', 'atom4']]
            impropers.atom1 = impropers.atom1.map(psf_types)
            impropers.atom2 = impropers.atom2.map(psf_types)
            impropers.atom3 = impropers.atom3.map(psf_types)
            impropers.atom4 = impropers.atom4.map(psf_types)
            impropers['atuple'] = list(zip(impropers.atom1, impropers.atom2, impropers.atom3, impropers.atom4))
            impropers.drop(columns=['atom1', 'atom2', 'atom3', 'atom4'], inplace=True)
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
            
            self.__init__(impropers)    
        	

#=====================================================================

class Region:
    '''
        Container for the Mixture.
    '''
    
    __DENSITY = {None:[0,0,1],"WATER": [0.9970479,298.15,18.01528], "Chloroform":[1.483,298.15,119.38], "Chloroform (4 atoms)": [1.483,298.15,119.38], "Acetone": [0.786,293.15,58.08], "THF": [0.886,293.15,72.11], "Argon": [5.537, 0.03049,39.948], "Methanol": [0.791,293.15,32.04], "Selected Molecule": [0.9970479,298.15,18.01528]}
    __DI = 0
    __TI = 1
    __MI = 2

    def __init__(self, cid=1):
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
            D = self.__DENSITY[solvent][self.__DI]
            MM = self.__DENSITY[solvent][self.__MI]
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
    def __init__(self, cid=1, maxsMins=None):
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



#===================================================================

#Aqui estaba las clases Df del MolecularTopolyData
        
#===================================================================

#Aqui estaba las las clases Coeffs del CharmmForceField()
 

#=================================================================== 

class LammpsData():
    ''' Holds LAMMPS data.
        See "Format of a data file" section in https://lammps.sandia.gov/doc/read_data.html
    '''

    def __init__(self,file=None, atoms=AtomsFull(), topology=MolecularTopologyData(), forceField=CharmmForceField(), region=Box()):

        '''
        self['Ellipsoids']  = EllipsoidsDF()
        selflines       = LinesDF()
        self.bodies      = BodiesDF()
        self.triangles   = TrianglesDF()
        self.bodies   = BodiesDF()
        '''
        
        #forceFieldSectio composi
        assert isinstance(forceField, CharmmForceField)
        self.forceField = forceField
        
        # atom-property sections
        assert isinstance(atoms, AtomsFull)
        self.atomprop = atoms
        
        # molecular topology sections
        assert isinstance(topology, MolecularTopologyData)
        self.topology = topology
        
        assert isinstance(region, Box)
        self.region = region
        
        if file:
            self.read(file)
        

    def read(self, filename):
        ''' Reads a LAMMS .data file.
        '''
        
        #Abrir el archivo para leer datos
        arch = open(filename, 'r')
        #serial = 1
        #ind = 0
        #num = []
        #tipo = []
        #caja = []
        
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
                if key == 'Masses': self.atomprop.masses.add(data)
                elif key == 'Pair Coeffs': self.forceField.pairCoeffs.add(data)
                elif key == 'Bond Coeffs': self.forceField.bondCoeffs.add(data)
                elif key == 'Angle Coeffs': self.forceField.angleCoeffs.add(data)                
                elif key == 'Dihedral Coeffs': self.topology.dihedralCoeffs.add(data)
                elif key == 'Improper Coeffs':self.forceField.improperCoeffs.add(data)
                elif key == 'Atoms': self.atomprop.atoms.add(data)                
                elif key == 'Velociteies': self.atomprop.velocities.add(data)   
                elif key == 'Bonds': self.topology.bonds.add(data)
                elif key == 'Angles': self.topology.angles.add(data)                
                elif key == 'Dihedrals': self.topology.dihedrals.add(data)
                elif key == 'Impropers':self.topology.impropers.add(data)                
              
        arch.close()                
           

    def loadNAMDdata(self, charmm):
        ''' loads data from NAMDdata object into self.'''
        #print("loadNAMDdata=",charmm.psf.dihedrals)
        
        #AtomPropertyData
        #if not isinstance(self.atomprop, AtomsFull): self.atomprop = AtomsFull()
        self.atomprop.setFromNAMD(charmm)
        
        # MolecularTopologyData
        #self.topology = MolecularTopologyData()
        self.topology.setFromNAMD(charmm)
       
        # CharmmForceField
        #self.forceField = CharmmForceField()
        self.forceField.setFromNAMD(charmm, self.atomprop,self.topology)
        
        # MolecularTopologyData
        #self.region = Box()
        self.region.setFromNAMD(charmm)
       
        '''
        self.pairCoeffs.setFromNAMD(charmm, self.masses)
        self.bondCoeffs.setFromNAMD(charmm, self.bonds)
        self.angleCoeffs.setFromNAMD(charmm, self.angles)
        self.dihedralCoeffs.setFromNAMD(charmm, self.dihedrals)
        self.improperCoeffs.setFromNAMD(charmm, self.impropers)
        '''

    def loadWolffia(self, wolffia, progressTitle=None):
        '''
        Converts a WolffiaState to a LammpsData self object.

        @param: mix a Wolffia Mixture object.
        '''
        from  granules.structure.NAMDdata import NAMDdata
        
        charmm = NAMDdata()
        print("LAMMPS.loadWolffia charmm.loadWolffia(wolffia)")
        charmm.loadWolffia(wolffia, progressTitle=progressTitle)
        print("LAMMPS.loadWolffia self.loadNAMDdata(charmm)")
        self.loadNAMDdata(charmm)
 
        return self


    def writeConf(self, filename, progressTitle=None):
        import sys
        from wolffialib.io.PrintBar import PrintBar
              
        if progressTitle != None:    
            progress = PrintBar(0, 13)
            progress.setLabelText(progressTitle+", writing conf")
            progress.setRange(0,13)
            progress.setValue(0)
        else:
            progress = None

        cfile = open(filename, "w")


        # Sección automática
        overview_section = '{0:>12}  atoms\n{1:>12}  bonds\n{2:>12}  angles\n{3:>12}  dihedrals\n{4:>12}  impropers\n\n{5:>12}  atom types\n' \
                           '{6:>12}  bond types\n{7:>12}  angle types\n{8:>12}  dihedral types\n{9:>12}  improper types\n\n' \
                           .format( \
                                len(self.atomprop.atoms), len(self.topology.bonds), len(self.topology.angles), len(self.topology.dihedrals), len(self.topology.impropers), len(self.atomprop.masses), \
                                len(self.forceField.bondCoeffs), len(self.forceField.angleCoeffs), len(self.forceField.dihedralCoeffs), len(self.forceField.improperCoeffs)) 

        try:
            maxminstr = [str(x) for x in self.region.maxsMins]
            box_section =  ' '.join(maxminstr[:2])  + ' xlo xhi\n' + \
                           ' '.join(maxminstr[2:4]) + ' ylo yhi\n' + \
                           ' '.join(maxminstr[4:])  + ' zlo zhi\n'
        except:
            # Improve and ask from where we get this data
            box_section =  ' ' + str(self.atomprop.atoms['x'].min()-2) + ' ' + str(self.atomprop.atoms['x'].max()+2) + ' xlo xhi\n' + \
                           ' ' + str(self.atomprop.atoms['y'].min()-2) + ' ' + str(self.atomprop.atoms['y'].max()+2) + ' ylo yhi\n' + \
                           ' ' + str(self.atomprop.atoms['z'].min()-2) + ' ' + str(self.atomprop.atoms['z'].max()+2) + ' zlo zhi\n'   
                       
        # Header
        cfile.write("LAMMPS Description\n\n")
        cfile.write(overview_section + box_section)
        if progress != None:    progress.setValue(1)


        # Molecule ID index
        if len(self.atomprop.atoms.moleculeTypes) > 0:
            cfile.write('\n# Molecule ID index\n#\n')
            for mID in self.atomprop.atoms.moleculeTypes:
                cfile.write("#{:10d}  {}\n".format(self.atomprop.atoms.moleculeTypes[mID],mID))
            cfile.write("\n")


        # Atom types index
        if len(self.atomprop.atoms.atomTypes) > 0:
            cfile.write('\n# Atom types index\n#\n')
            for mID in self.atomprop.atoms.atomTypes:
                cfile.write("#{:10d}  {}\n".format(mID, self.atomprop.atoms.atomTypes[mID]))
            cfile.write("\n")

        #Masses
        if len(self.atomprop.masses) > 0:
            cfile.write('\nMasses\n\n')
            cfile.write(self.atomprop.masses.to_string(index=False, columns=self.atomprop.masses.columns, header=False))
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No atom mass values to write. Simulation unlikely to run with this file.")
        if progress != None:    progress.setValue(2)


        #PairIJ Coeffs
        if len(self.forceField.pairCoeffs) > 0:
            cfile.write('\nPairIJ Coeffs\n\n')
            #print("Pair Coeffs:", self.pairCoeffs.columns)
            #[cfile.write('{:>3d}{:>12}{:>12}{:>12}{:>12}\n'.format(row['ID'], row['Charge'], row['Energy'],
            # row['Charge'], row['Energy'])) for index, row in self['Pair Coeffs'].iterrows()]
            cfile.write(self.forceField.pairCoeffs.to_string(index=False, columns=self.forceField.pairCoeffs.columns, header=False))        
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No Pair coefficients to write.\n")
        if progress != None:    progress.setValue(3)


        #Bond Coeffs
        if len(self.forceField.bondCoeffs) > 0:
            cfile.write('\nBond Coeffs\n\n')
            cfile.write(self.forceField.bondCoeffs.to_string(index=False, columns=self.forceField.bondCoeffs.columns, header=False))        
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No bond coefficients to write.\n")
        if progress != None:    progress.setValue(4)


        #Angle Coeffs
        if len(self.forceField.angleCoeffs) > 0:
            cfile.write('\nAngle Coeffs\n\n')
            cfile.write(self.forceField.angleCoeffs.to_string(index=False, columns=self.forceField.angleCoeffs.columns, header=False))        
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No angle coefficients to write.\n")
        if progress != None:    progress.setValue(5)


        #Dihedral Coeffs
        if len(self.forceField.dihedralCoeffs) > 0:
            cfile.write('\nDihedral Coeffs\n\n')
            cfile.write(self.forceField.dihedralCoeffs.to_string(index=False, columns=self.forceField.dihedralCoeffs.columns, header=False))        
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No dihedral coefficients to write.\n")
        if progress != None:    progress.setValue(6)


        #Improper Coeffs
        if len(self.forceField.improperCoeffs) > 0:
            cfile.write('\nImproper Coeffs\n\n') 
            cfile.write(self.forceField.improperCoeffs.to_string(index=False, columns=self.forceField.improperCoeffs.columns, header=False))        
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No improper coefficients to write.\n")
        if progress != None: progress.setValue(7)


        #Atoms
        cfile.write('\nAtoms\n\n') 
        #cfile.write(self.atomprop.atoms.to_string(index=False, columns=self.atomprop.atoms.columns, header=False))
        cfile.write(self.atomprop.atoms.to_string(index=False, columns=self.atomprop.atoms.columns, header=False))
        cfile.write("\n")
        if progress != None: progress.setValue(8)


        #Velocities
        if len(self.atomprop.velocities) > 0:
            cfile.write('\nVelocities\n\n')
            cfile.write(self.atomprop.velocities.to_string(index=False, columns=self.atomprop.velocities.columns, header=False))        
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No velocities to write.\n")
        if progress != None: progress.setValue(9)


        #Bonds
        if len(self.topology.bonds) > 0:
            cfile.write('\nBonds\n\n')
            cfile.write(self.topology.bonds.to_string(index=False, columns=self.topology.bonds.columns, header=False))
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No bonds to write.\n")
        if progress != None: progress.setValue(10)

        #Angles
        if len(self.topology.angles) > 0:
            cfile.write('\nAngles\n\n') 
            cfile.write(self.topology.angles.to_string(index=False, columns=self.topology.angles.columns, header=False))
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No angles to write.\n")
        if progress != None: progress.setValue(11)


        #Dihedrals
        if len(self.topology.dihedrals) > 0:
            cfile.write('\nDihedrals\n\n') 
            cfile.write(self.topology.dihedrals.to_string(index=False, columns=self.topology.dihedrals.columns, header=False))
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No dihedrals to write.\n")
        if progress != None: progress.setValue(12)


        #Impropers
        if len(self.topology.impropers) > 0:
            cfile.write('\nImpropers\n\n') 
            cfile.write(self.topology.impropers.to_string(index=False, columns=self.topology.impropers.columns, header=False))
            cfile.write("\n")
        else:
            sys.stderr.write("WARNING: No impropers to write.\n")
        if progress != None:
            progress.setValue(13)
            progress.close()
    
    
    def Force(self):
        '''Hace una llamada a la funcion Force() de la clase forceField() 
            para poder printiar los datos.'''
        
        print("LammpsData.Force()")
        self.forceField.Force(self.atomprop,self.topology)

    def append(self,other):
        '''Une dos objetos de LammpsData, sus dataframes individuales'''
        # OH YEAH
        self.forceField.angleCoeffs = self.forceField.angleCoeffs.append(other.forceField.angleCoeffs)
        self.forceField.bondCoeffs = self.forceField.bondCoeffs.append(other.forceField.bondCoeffs)
        self.forceField.dihedralCoeffs = self.forceField.dihedralCoeffs.append(other.forceField.dihedralCoeffs)
        self.forceField.improperCoeffs = self.forceField.improperCoeffs.append(other.forceField.improperCoeffs)
        self.forceField.pairCoeffs = self.forceField.pairCoeffs.append(other.forceField.pairCoeffs)
        #Oh YEaSH
        self.atomprop.atoms = self.atomprop.atoms.append(other.atomprop.atoms)
        self.atomprop.velocities = self.atomprop.velocities.append(other.atomprop.velocities)
        self.atomprop.masses = self.atomprop.masses.append(other.atomprop.masses)
        #yeyeyeye
        self.topology.angles = self.topology.angles.append(other.topology.angles)
        self.topology.bonds = self.topology.bonds.append(other.topology.bonds)
        self.topology.dihedrals = self.topology.dihedrals.append(other.topology.dihedrals)
        self.topology.impropers = self.topologtopologiaia.impropers.append(other.topology.impropers)


    '''        
    def copy(self):#modifica
        ld = LammpsData()
        # OH YEAH
        ld.forceField.angleCoeffs = self.forceField.angleCoeffs.copy()
        ld.forceField.bondCoeffs = self.forceField.bondCoeffs.copy()
        ld.forceField.dihedralCoeffs = self.forceField.dihedralCoeffs.copy()
        ld.forceField.improperCoeffs = self.forceField.improperCoeffs.copy()
        ld.forceField.pairCoeffs = self.forceField.pairCoeffs.copy()
        # OH YEAHHH
        ld.atomprop.atoms = self.atomprop.atoms.copy()
        ld.atomprop.velocities = self.atomprop.velocities.copy()
        ld.atomprop.masses = self.atomprop.masses.copy()
        # OH YEAHHH
        ld.topology.angles = self.topology.angles.copy()
        ld.topology.bonds = self.topology.bonds.copy()
        ld.topology.dihedrals = self.topology.dihedrals.copy()
        ld.topology.impropers = self.topology.impropers.copy()
        
        return ld
    '''
    
    def selectAtom(self,atomNumber):
        '''Funcion que elimina un tipo de atomo deseado del dataframe'''
        
        # AtompropertyData 
        #atoms
        self.atomprop.atoms.drop(self.atomprop.atoms[self.atomprop.atoms.aType == atomNumber].index, inplace=True)#elimina el atomo indicado
        self.atomprop.atoms.reset_index(inplace=True)#reset el index, crea una copia del index
        self.atomprop.atoms.drop(['index','aID'],axis = 1, inplace=True)#elimina la copia y la columna afectada
        self.atomprop.atoms.insert(0, "aID", [*range(1,len(self.atomprop.atoms)+1)], True)#inserta columna nueva en la afectada("aID")
        #velocity
        self.atomprop.velocities.drop(self.atomprop.velocities.tail(166-len(self.atomprop.atoms)).index,inplace = True) #elimina velocidades dependiendo a los atomos eliminados
        #masses
         #chequea
         
        #TopologiaData
        #bonds
        #Cambio para que busque en los Atomos enves de Type
        self.topology.bonds = self.topology.bonds.ix[self.topology.bonds['atom1'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topology.bonds = self.topology.bonds.ix[self.topology.bonds['atom2'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topology.bonds  = self.topology.bonds.reset_index(drop=True)#resetea el indice y elimina la copia
        self.topology.bonds['bID'] = [i for i in range(1,len(self.topology.bonds)+1 )]#resetea el indice en la columna 'bType'
        #angles
        self.topology.angles = self.topology.angles.ix[self.topology.angles['atom3'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topology.bonds = self.topology.bonds.ix[self.topology.bonds['atom1'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topology.bonds = self.topology.bonds.ix[self.topology.bonds['atom2'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topology.angles  = self.topology.angles.reset_index(drop=True)#resetea el indice y elimina la copia
        self.topology.angles['anID'] = [i for i in range(1,len(self.topology.angles)+1 )]#resetea el indice en la columna 'anType'
        #improper######modifica desde aqui
        self.topology.impropers = self.topology.impropers.ix[self.topology.impropers['iType'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topology.impropers  = self.topology.impropers.reset_index(drop=True)#resetea el indice y elimina la copia
        self.topology.impropers['iID'] = [i for i in range(1,len(self.topology.impropers)+1 )]#resetea el indice en la columna 'iType'
        #dihedral
        self.topology.dihedrals = self.topology.dihedrals.ix[self.topology.dihedrals['dType'] != atomNumber ]#devuelve los rows que no contengan el atomNumber
        self.topology.dihedrals  = self.topology.dihedrals.reset_index(drop=True)#resetea el indice y elimina la copia
        self.topology.dihedrals['dID'] = [i for i in range(1,len(self.topology.dihedrals)+1 )]#resetea el indice en la columna 'dType'
        
        
    def bondLength(self,bondsID=set()):
        ''' Calculates the length of the bond between two atoms '''

        if set(bondsID) == set(): 
            bondsID = set(self.topology.bonds.bID)
            
        coordinates = self.atomprop.atoms[['aID', 'x', 'y', 'z']].set_index('aID')
        select = self.topology.bonds.loc[self.topology.bonds['bID'].isin(bondsID)][['bID', 'bType','atom1', 'atom2']]
        #print("bondLength coordinates {}\n{}".format(type(coordinates),coordinates))
        #print("bondLength select {}\n{}".format(type(select),select))
        #print("bondLength select.join(coordinates, on='atom1') {}\n{}".format(type(select.join(coordinates, on='atom1')),select.join(coordinates, on='atom1')))
        a1 = select.join(coordinates, on='atom1')[['x','y','z']]
        a2 = select.join(coordinates, on='atom2')[['x','y','z']]
        #reset_index para iterar por los valores al combinar
        select['dist'] = (np.sqrt(((a1-a2)**2).sum(axis=1)))#.reset_index(drop=True)
        return select[['bID', 'bType','dist']].set_index('bID')

    def meanBondLengthByType(self):
        '''Calculates the average distance of the Bond types'''
        bLengths = self.bondLength()
        bLengths['bType'] = self.topology.bonds.bType
        return bLengths.groupby('bType').mean()

    def describeBonds(self):
        ''' return the descriptive statistics by bond type.
        '''
        return self.bondLength().groupby('bType').describe()

    def angleLength(self,atomIds=set()):
        '''Calculates the length of the angles between 3 atoms'''
        if set(atomIds) == set(): 
            atomIds = set(self.topology.angles.anID)        
        coordinates = self.atomprop.atoms[['aID', 'x', 'y', 'z']].set_index('aID')
        selection = self.topology.angles.loc[self.topology.angles['anID'].isin(atomIds)].copy() 
        print("angleLength", type(selection))
        a1 = selection.join(coordinates, on='atom1')[['x','y','z']]
        a2 = selection.join(coordinates, on='atom2')[['x','y','z']]
        a3 = selection.join(coordinates, on='atom3')[['x','y','z']]
        #print(self.topology.angles.loc[self.topology.angles['atom3']==3])
        #print("angleLength: " , a1,"\n",a2,"\n",a3)
        
        #vectores de ambos lados 
        v1 = a1 - a2
        v2 = a3 - a2
        #print('#######calculos',v1,'\n',v2)
        #sacamos el dot product(multiplicar vectores)
        d = np.einsum('ij,ij->i', v1, v2) 
        #print("anglelentgth: v1 * v2", d)
        #calcula la magnitud de los lados
        absolV1 = np.sqrt(((v1)**2).sum(axis=1))
        absolV2 = np.sqrt(((v2)**2).sum(axis=1))
        
        #print("anglelength abs:", absolV1, absolV2)
        #Calculamos angulo
        product = d/(absolV1*absolV2)
        #print("anglelength return:", np.arccos(product))
        selection['dist'] = np.arccos(product).reset_index(drop=True)
        return selection[['anType','dist']]
    
    def meanAngleLengthByType(self):
        '''Calculates the average distance of the angles types'''
        anLengths = self.angleLength()
        anLengths['anType'] = self.topology.angles.anType
        return anLengths.groupby('anType').mean()
    
    def describeAngle(self):
        ''' return the descriptive statistics by angle type.
        '''
        return self.angleLength().groupby('anType').describe()

    
    ''' Esto crea información redundante con riesgo de que se vuelva inconsistente
    def combineLength(self,bondTable,angleTable):
        #Stores the lenght of bond and angle in single panda table
        
        table = pd.DataFrame({'Column_0':[bondTable.loc[0], angleTable.loc[0]]})
        #print("aqui#2")
        for index in range(1,len(bondTable)):
            table['Column_'+str(index)] = [bondTable.loc[index],angleTable.loc[index]]
            #table = table.append({'IDs': bondTable.loc[index]}, ignore_index = True).reset_index(drop=True)
            #table = table.append({'IDs':angleTable.loc[index]}, ignore_index = True).reset_index(drop=True)

            #table.append(angleTable.loc[index]).reset_index(drop=True)
        return table 
       
               
        #return bondTable.copy().append(angleTable).reset_index(drop=True)
     '''
_AVOGRADRO_CONSTANT_ = 6.02214129e+23



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
