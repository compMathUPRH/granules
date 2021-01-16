# -*- coding: utf-8 -*-
"""-------------------------------------------------------------------------
  NAMDdata.py
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


class PDB(pd.DataFrame):
    ''' Pandas DataFrame that stores data defined by the PDB format 
        (only the ATOMS section).
    '''

    def __init__(self,data=None, dtype=None, copy=False):
        super(PDB, self).__init__(data=data, 
                                  copy=copy,
                                  columns=['RecName','ID','Name','AltLoc','ResName',
                                           'ChainID','ResSeq','iCode','x','y','z',
                                           'Occupancy','TempFactor','Element','Charge']
        )


    def select(self, selection, reverse=False):
        '''
        try:
            col = 'Name'
            self.__init__(data=self[self[col].isin(selection[col])])
        except: pass
        '''
        for col in selection:
            if col in self.columns:
                if reverse:
                    self.__init__(data=self[~self[col].isin(selection[col])])
                else:
                    self.__init__(data=self[self[col].isin(selection[col])])
            else:
                print("WARNING: Selection index \"{}\" was not recognized and it was ignored. Valid options are {}".format(col, self.columns))


    def readFile(self, filename):
        ''' reads PDB file and appends to self.
            Follows specifications of the atoms section in 
            http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

            Parameter
            ----------
            filename : str
                name of file
        '''

        #print("PDB.readFile(",filename,")")
        
        #Abrir el archivo para leer datos
        arch = open(filename, 'r')
        data = []
        
        ranges = [6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6,10,2,2]
        voidPositions = [0,6,11,20,27,66]
        
        atomID = 1
        for linea in arch:
            
            #Solo queremos la informacion de átomos
            if linea[:4] == 'ATOM' or linea[:4] == 'HETA': 
                info = ['ATOM',atomID]
                start = 0
                stop = 0
                        
                for i in ranges:
                    stop += i
                    if start not in voidPositions:
                        l = linea[start:stop].strip()
                        #print('\nEn posicion',start+1,':',stop,'se encuentra ->>',l)
                        if l == '' or l == '<0>': info.append(np.nan)
                        else: info.append(l)
                    start += i

                #Cada vez anadir nueva informacion a los datos
                #print(info)
                data.append(info)
                atomID += 1
                       
        arch.close() 
        newTable = PDB(data=data)  
        #print(newTable)

        # set column types and append to existing table
        super().__init__(data=self.append(newTable, ignore_index=True).astype(
                {'ID'        :int,
                 'ResSeq'    :int,
                 'x'         :float,
                 'y'         :float,
                 'z'         :float,
                 'Occupancy' :float,
                 'TempFactor':float  
                 #'Charge'    :float 
                })
        )

        #print("PDB.readFile(",filename,") ... END")
            
class PBC():
    ''' Stores periodic boundary conditions as 3 cell vectors and the coordinates
        of the origin.
        
        Parameters;
        cellv1, cellv, cellv3 and cello: np.arrays or lists with 3 float values
        : 
        
    '''

    def __init__(self,cellv1=None, cellv2=None, cellv3=None, cello=None, 
                 wrapWater=False, wrapNearest=False):
        # validate parameters
        A = cellv1==None or cellv1==None or cellv1==None
        B = cellv1==None and cellv1==None and cellv1==None
        assert(not A or (A and B) or (A and not B)) # A implies B
        
        self.cellBasisVector1 = cellv1
        self.cellBasisVector2 = cellv2
        self.cellBasisVector3 = cellv3
        self.cellOrigin       = cello
        self.wrapWater        = wrapWater
        self.wrapNearest      = wrapNearest

    def readFile(self, filename):
        ''' reads cell basis vectors and cell origin from XSC file. Ignores the rest. 

            Parameter
            ----------
            filename : str
                name of file
        '''
        import numpy as np
        
        #print("PBC.readFile(",filename,")")
        xsc_file = open(filename, "r")
        for line in xsc_file:
            #print("PBC.readFile(",filename,") ... line = ", line)
            if line.strip()[0] != '#':
                record = line.strip().split(' ')
                self.cellBasisVector1 = np.array(record[1:4],   dtype='float')
                self.cellBasisVector2 = np.array(record[4:7],   dtype='float')
                self.cellBasisVector3 = np.array(record[7:10],  dtype='float')
                self.cellOrigin       = np.array(record[10:13], dtype='float')
                #print("PBC.readFile(",filename,") ... record = ", record)
                #print("PBC.readFile(",filename,") ... self.cellBasisVector1 = ", self.cellBasisVector1)
                #print("PBC.readFile(",filename,") ... self.cellBasisVector2 = ", self.cellBasisVector2)
                #print("PBC.readFile(",filename,") ... self.cellBasisVector3 = ", self.cellBasisVector3)
                #print("PBC.readFile(",filename,") ... self.cellOrigin = ", self.cellOrigin)
        #print("PBC.readFile(",filename,") ... END")
       
#==============================================================================
class PSF:
    ''' Pandas DataFrames that store data defined in the PSF file format specificaton.
        See specification in 
        http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html
    '''

    def __init__(self):
        self.atoms = PSF.ATOM()
        self.bonds = PSF.BOND()
        self.angles = PSF.THETA()
        self.dihedrals = PSF.PHI()
        self.impropers = PSF.IMPHI()
        self.cross_terms = PSF.CRTERM()

    def readFile(self, filename):
        ''' Reads data for all the sections in the PSF intto self.

            Parameters:
            -------------------
            filename : str
                psf file name
        '''
        #print("PSF.readFile(",filename,")")
        
        self.atoms.readSection(filename)
        self.bonds.readSection(filename)
        self.angles.readSection(filename)
        self.dihedrals.readSection(filename)
        self.impropers.readSection(filename)
        self.cross_terms.readSection(filename)
        
        #print("PSF.readFile(",filename,") ... END")

    def select(self, atoms):
        for table in [self.atoms, self.bonds, self.angles, self.dihedrals, self.impropers, self.cross_terms]:
            table.select(atoms)

    @staticmethod
    def readSection(filename, section, tupleLength, itemsPerLine):
        ''' reads the section of the PSF file specified in the parameter
             'section'.

            Parameter
            ----------
            filename : str
                name of file

            section : str
                one of "ATOM", "BOND", ...

            itemsPerLine : int
                hoy many bonds, angles, ... are in each line

            Returns
                lists of lists of values
        '''
        
        #Abrir el archivo para leer datos
        arch = open(filename, 'r')
        data = []        

        # find desired section
        for linea in arch:
            if '!N'+section in linea: break

        #Obtener la cantidad de elementos, cae a 0 si no encontró la sección 
        try:
            cantidad = int(linea[:9].strip())
        except:
            cantidad = 0    
        
        if cantidad % itemsPerLine == 0: cantLineas = cantidad // itemsPerLine
        else: cantLineas = cantidad // itemsPerLine + 1
        #print("cantLineas {0}: {1}".format(section,cantLineas))

        for i in range(cantLineas):
            linea = next(arch)
            #Recoger la informacion mediante su posicion
            info = linea.split()
            for j in range(0,len(info), tupleLength):
                #Cada vez anadir nueva informacion a los datos
                data.append(info[j:j+tupleLength])  
            #print(">", info)

        arch.close()
        #print(data)
        return data
    

    #---------------- inner class of PSF --------------------------------------
    class ATOM(pd.DataFrame):
        ''' ATOM section of the PSF file format specification.'''

        def __init__(self,data=None, dtype=None, copy=False):
            super(PSF.ATOM, self).__init__(data=data, copy=copy, columns=[
                        'ID','RecName','ChainID', 'ResName', 'Name', 
                        'Type', 'Charge', 'Mass', 'Unused'])

        def readSection(self, filename):
            ''' reads the section of the ATOM section of PSF file specified in the parameter
                 'filename'.
    
                Parameter
                ----------
                filename : str
                    name of file
            '''
            newTable = PSF.ATOM(data=PSF.readSection(filename, "ATOM", 9, 1))  

            # set column types and append to existing table
            super().__init__(data=self.append(newTable, ignore_index=True).astype({
                     'ID'     :int,
                     'Charge' :float,
                     'Mass'   :float
                    }))

        def select(self, atoms):
            self.__init__(data=self[self['ID'].isin(atoms)])



    #---------------- inner class of PSF --------------------------------------
    class BOND(pd.DataFrame):
        ''' BOND section of the PSF file format specification.'''

        def __init__(self,data=None, dtype=None, copy=False):
            super(PSF.BOND, self).__init__(data=data, copy=copy, columns=[
                    'atom1','atom2'])

        def readSection(self, filename):
            ''' reads the section of the BOND section of PSF file specified in the parameter
                 'filename'.
    
                Parameter
                ----------
                filename : str
                    name of file
            '''
            data=PSF.readSection(filename, "BOND", 2, 4)
            newTable = pd.DataFrame(data).dropna()
            pair1 = newTable[[0,1]]

            pair1.columns=['atom1','atom2']

            super().__init__(data=self.append(pair1, ignore_index=True).astype({
                     'atom1'     :int,
                     'atom2' :int
                    }))

        def select(self, atoms):
            self.__init__(data=self[self['atom1'].isin(atoms) & (self['atom2'].isin(atoms))])

    #---------------- inner class of PSF --------------------------------------
    class THETA(pd.DataFrame):
        ''' THETA section of the PSF file format specification.'''

        def __init__(self,data=None, dtype=None, copy=False):
            super(PSF.THETA, self).__init__(data=data, copy=copy, columns=[
                    'atom1','atom2','atom3'])

        def readSection(self, filename):
            ''' reads the section of the THETA section of PSF file specified in the parameter
                 'filename'.
    
                Parameter
                ----------
                filename : str
                    name of file
            '''
            data=PSF.readSection(filename, "THETA", 3,3)
            newTable = pd.DataFrame(data).dropna()
            pair1 = newTable[[0,1,2]]
            pair1.columns=['atom1','atom2','atom3']

            super().__init__(data=self.append(pair1, ignore_index=True).astype({
                     'atom1'     :int,
                     'atom2'     :int,
                     'atom3' :int
                    }))

        def select(self, atoms):
            self.__init__(data=self[self['atom1'].isin(atoms) & (self['atom2'].isin(atoms)) & (self['atom3'].isin(atoms))])

    #---------------- inner class of PSF --------------------------------------
    class PHI(pd.DataFrame):
        ''' PHI section of the PSF file format specification.'''

        def __init__(self,data=None, dtype=None, copy=False):
            super(PSF.PHI, self).__init__(data=data, copy=copy, columns=[
                    'atom1','atom2','atom3','atom4'])

        def readSection(self, filename):
            ''' reads the section of the PHI section of PSF file specified in the parameter
                 'filename'.
    
                Parameter
                ----------
                filename : str
                    name of file
            '''
            data=PSF.readSection(filename, "PHI", 4,2)
            if len(data) == 0:
                super().__init__(data=self.astype({
                     'atom1'     :int,
                     'atom2'     :int,
                     'atom3'     :int,
                     'atom4' :int
                    }))
            else:
                newTable = pd.DataFrame(data).dropna()
                pair1 = newTable[[0,1,2,3]]
                pair1.columns=['atom1','atom2','atom3','atom4']

                super().__init__(data=self.append(pair1, ignore_index=True).astype({
		                 'atom1'     :int,
		                 'atom2'     :int,
		                 'atom3'     :int,
		                 'atom4' :int
		                }))


        def select(self, atoms):
            self.__init__(data=self[self['atom1'].isin(atoms) & (self['atom2'].isin(atoms)) & (self['atom3'].isin(atoms)) & (self['atom4'].isin(atoms))])


    #---------------- inner class of PSF --------------------------------------
    class IMPHI(pd.DataFrame):
        ''' IMPHI section of the PSF file format specification.'''

        def __init__(self,data=None, dtype=None, copy=False):
            super(PSF.IMPHI, self).__init__(data=data, copy=copy, columns=[
                    'atom1','atom2','atom3','atom4'])

        def readSection(self, filename):
            ''' reads the section of the IMPHI section of PSF file specified in the parameter
                 'filename'.
    
                Parameter
                ----------
                filename : str
                    name of file
            '''
            data=PSF.readSection(filename, "IMPHI", 4,2)
            if len(data) == 0:
                super().__init__(data=self.astype({
                     'atom1'     :int,
                     'atom2'     :int,
                     'atom3'     :int,
                     'atom4' :int
                    }))
            else:
                newTable = pd.DataFrame(data).dropna()
                pair1 = newTable[[0,1,2,3]]
                pair1.columns=['atom1','atom2','atom3','atom4']

                super().__init__(data=self.append(pair1, ignore_index=True).astype({
                         'atom1'     :int,
                         'atom2'     :int,
                         'atom3'     :int,
                         'atom4' :int
                        }))

        def select(self, atoms):
            self.__init__(data=self[self['atom1'].isin(atoms) & (self['atom2'].isin(atoms)) & (self['atom3'].isin(atoms)) & (self['atom4'].isin(atoms))])


    #---------------- inner class of PSF --------------------------------------
    class CRTERM(pd.DataFrame):
        ''' CRTERM section of the PSF file format specification.'''

        def __init__(self,data=None, dtype=None, copy=False):
            super(PSF.CRTERM, self).__init__(data=data, copy=copy, columns=[
                    'atom1','atom2','atom3','atom4'])

        def readSection(self, filename):
            ''' reads the section of the CRTERM section of PSF file specified in the parameter
                 'filename'.
    
                Parameter
                ----------
                filename : str
                    name of file
            '''
            data=PSF.readSection(filename, "CRTERM", 4,2)
            if len(data) == 0:
                super().__init__(data=self.astype({
                     'atom1'     :int,
                     'atom2'     :int,
                     'atom3'     :int,
                     'atom4' :int
                    }))
            else:
                newTable = pd.DataFrame(data).dropna()
                pair1 = newTable[[0,1,2,3]]
                pair1.columns=['atom1','atom2','atom3','atom4']
                super().__init__(data=self.append(pair1, ignore_index=True))

        def select(self, atoms):
            self.__init__(data=self[self['atom1'].isin(atoms) & (self['atom2'].isin(atoms)) & (self['atom3'].isin(atoms)) & (self['atom4'].isin(atoms))])

#==============================================================================
class PRM:
    ''' NAMDdata force field from a PRM file (CHARMM force field parameters).
        See format in https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node25.html
    '''
    #SECTIONS = ["BONDS", "ANGLES", "DIHEDRALS", "IMPROPER", "CMAP", "NONBONDED", "END", "HBOND"]
    SECTIONS = ["BOND", "ANGL", "DIHE", "IMPR", "CMAP", "NONB", "END", "HBON"]

    def __init__(self):
        self.nonbonded = PRM.NONBONDED()
        self.bonds = PRM.BONDS()
        self.angles = PRM.ANGLES()
        self.dihedrals = PRM.DIHEDRALS()
        self.impropers = PRM.IMPROPER()


    def readFile(self, filename):
        ''' Reads data for all the sections in the PRM .

            Parameters:
            -------------------
            filename : str
                prm file name
        '''

        #print("PRM.readFile(",filename,")")
        self.bonds.readSection(filename)
        self.angles.readSection(filename)
        self.dihedrals.readSection(filename)
        self.impropers.readSection(filename)
        self.nonbonded.readSection(filename)
        #print("PRM.readFile(",filename,") ... END")

   
    @staticmethod
    def readSection(filename, section):
        ''' reads the section of the PRM file specified in the parameter
             'section'.

            Parameter
            ----------
            filename : str
                name of file

            section : str
                a string in PRM.SECTIONS

            Returns
                lists of lists of string values
        '''
        
        #Abrir el archivo para leer datos
        arch = open(filename, 'r')
        data = []       

        # find desired section
        while True:  # search for multiple sections, exits when it hits the end of file
            found = False
            for linea in arch:
                linea = linea.split("!")[0].strip()
                if len(linea) == 0 or linea[0] == '*': continue  # skip comment
                if section[:4] == linea[:4]: 
                    found = True
                    break
            #print("PRM linea encontrada: ", linea)
            if not found: 
                arch.close()
                return data

            if linea[-1] == '-': # skips continuation line
                linea = arch.readline().split("!")[0].strip()

            for linea in arch:
                linea = linea.split("!")[0].strip()
                if len(linea) == 0 or linea[0] == '*': continue  # skip comment
                #print(">"+linea+'|',linea in PRM.SECTIONS)
                if linea.split(' ')[0][:4] in PRM.SECTIONS: break  # start of next section
                if linea != "": # skip empty or comment_only lines
                    info = linea.split()
                    #print(len(info), info)
                    data.append(info)  

        arch.close()
        return data


    #---------------- inner class of PRM --------------------------------------
    class NONBONDED(pd.DataFrame):
        ''' Stores information from the NONBONDED section of the 
            NAMDdata force field from a PRM file (CHARMM force field parameters).
            See format in https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node25.html
        '''
        def __init__(self,data=None, dtype=None, copy=False):
            super(PRM.NONBONDED, self).__init__(data=data, copy=copy, columns=[
                'Type','epsilon','Rmin2','epsilon1_4','Rmin2_1_4'])

        def getCoeffs(self):
            # retrieve info from PRM
            prmFF = self.copy()
            prmFF.set_index('Type', inplace=True)
            return prmFF


        def readSection(self, filename):
            ''' reads the NONBONDED section of the PRM file specified in the parameter
                 'section'.
    
                Parameter
                ----------
                filename : str
                    name of file
            '''

            data = list()

            d = PRM.readSection(filename, "NONBONDED")

            for list_ in d:
                #print(list_)
                if len(list_) < 7:
                    data.append([list_[0], -float(list_[2]), list_[3], -float(list_[2]), list_[3]])
                else:
                    data.append([list_[0], -float(list_[2]), list_[3], -float(list_[5]), list_[6]])

            newTable = PRM.NONBONDED(data=data)  

            # modificar newTable para que todas las listas tengan 4 números
            #print(newTable)

            super().__init__(data=self.append(newTable, ignore_index=True).astype({
                     'epsilon' :float,
                     'Rmin2' :float,
                     'epsilon1_4' :float,
                     'Rmin2_1_4' :float
                    }))




    #---------------- inner class of PRM --------------------------------------
    class BONDS(pd.DataFrame):
        def __init__(self,data=None, dtype=None, copy=False):
            super(PRM.BONDS, self).__init__(data=data, copy=copy, columns=[
                'Type1','Type2','Kb','b0'])

        def getCoeffs(self):
            # retrieve info from PRM
            prmFF = self.copy()
            prmFF['atuple'] = list(zip(prmFF.Type1, prmFF.Type2))
            prmFF.drop(columns=['Type1', 'Type2'], inplace=True)
            prmFF2 = self.copy()
            prmFF2['atuple'] = list(zip(prmFF2.Type2, prmFF2.Type1))  # flipped tupples
            prmFF2.drop(columns=['Type1', 'Type2'], inplace=True)
            prmFF = prmFF.append(prmFF2)
            prmFF.set_index('atuple', inplace=True)

            return prmFF


        def readSection(self, filename):
            ''' reads the BONDS section of the PRM file specified in the parameter
                 'section'.
    
                Parameter
                ----------
                filename : str
                    name of file
            '''
            newTable = PRM.BONDS(data=PRM.readSection(filename, "BONDS"))  
            super().__init__(data=self.append(newTable, ignore_index=True).astype({
                     'Kb' :float,
                     'b0'   :float
                    }))


    #---------------- inner class of PRM --------------------------------------
    class ANGLES(pd.DataFrame):
        def __init__(self,data=None, dtype=None, copy=False):
            super(PRM.ANGLES, self).__init__(data=data, copy=copy, columns=[
                'Type1','Type2','Type3','Ktheta','Theta0','Kub','S0'])

        def getCoeffs(self):
            # retrieve info from PRM
            prmFF = self.copy()
            prmFF['atuple'] = list(zip(prmFF.Type1, prmFF.Type2, prmFF.Type3))
            prmFF.drop(columns=['Type1', 'Type2', 'Type3'], inplace=True)
            prmFF2 = self.copy()
            prmFF2['atuple'] = list(zip(prmFF2.Type3, prmFF2.Type2, prmFF2.Type1))  # flipped tupples
            prmFF2.drop(columns=['Type1', 'Type2', 'Type3'], inplace=True)
            prmFF = prmFF.append(prmFF2)
            prmFF.set_index('atuple', inplace=True)

            return prmFF

        def readSection(self, filename):
            ''' reads the ANGLES section of the PRM file specified in the parameter
                 'section'.
    
                Parameter
                ----------
                filename : str
                    name of file
            '''
            data=PRM.readSection(filename, "ANGLES")

            if data != []:
                for row in data:
                    if len(row) == 5:  # add two columns to data
                        row +=  [np.nan, np.nan]
                        
            newTable = PRM.ANGLES(data=data)  
            super().__init__(data=self.append(newTable, ignore_index=True).astype({
                     'Ktheta' :float,
                     'Theta0'   :float,
                     'Kub' : float,
                     'S0' : float
                    }))


    #---------------- inner class of PRM --------------------------------------
    class DIHEDRALS(pd.DataFrame):
        def __init__(self,data=None, dtype=None, copy=False):
            super(PRM.DIHEDRALS, self).__init__(data=data, copy=copy, columns=[
                'Type1','Type2','Type3','Type4','Kchi',
                'n','delta'])

        def getCoeffs(self):
            # retrieve info from PRM
            prmFF = self.copy()
            prmFF['atuple'] = list(zip(prmFF.Type1, prmFF.Type2, prmFF.Type3, prmFF.Type4))
            prmFF.drop(columns=['Type1', 'Type2', 'Type3', 'Type4'], inplace=True)
            prmFF2 = self.copy()
            prmFF2['atuple'] = list(zip(prmFF2.Type4, prmFF2.Type3, prmFF2.Type2, prmFF2.Type1))  # flipped tupples
            prmFF2.drop(columns=['Type1', 'Type2', 'Type3', 'Type4'], inplace=True)
            prmFF = prmFF.append(prmFF2)
            prmFF.set_index('atuple', inplace=True)

            return prmFF

        def readSection(self, filename):
            ''' reads the DIHEDRALS section of the PRM file specified in the parameter
                 'section'.
    
                Parameter
                ----------
                filename : str
                    name of file
            '''
            newTable = PRM.DIHEDRALS(data=PRM.readSection(filename, "DIHEDRALS"))  
            super().__init__(data=self.append(newTable, ignore_index=True).astype({
                     'Kchi'  :float,
                     'n'     :int,
                     'delta' : float
                    }))

    #---------------- inner class of PRM --------------------------------------
    class IMPROPER(pd.DataFrame):
        def __init__(self,data=None, dtype=None, copy=False):
            super(PRM.IMPROPER, self).__init__(data=data, copy=copy, columns=[
                'Type1','Type2','Type3','Type4','Kpsi','unused','psi0'])

        def getCoeffs(self):
            # retrieve info from PRM
            prmFF = self.copy()
            prmFF['atuple'] = list(zip(prmFF.Type1, prmFF.Type2, prmFF.Type3, prmFF.Type4))
            prmFF.drop(columns=['Type1', 'Type2', 'Type3', 'Type4'], inplace=True)
            prmFF2 = self.copy()
            prmFF2['atuple'] = list(zip(prmFF2.Type4, prmFF2.Type3, prmFF2.Type2, prmFF2.Type1))  # flipped tupples
            prmFF2.drop(columns=['Type1', 'Type2', 'Type3', 'Type4'], inplace=True)
            prmFF = prmFF.append(prmFF2)
            prmFF.set_index('atuple', inplace=True)

            return prmFF

        def readSection(self, filename):
            ''' reads the IMPROPER section of the PRM file specified in the parameter
                 'section'.
    
                Parameter
                ----------
                filename : str
                    name of file
            '''
            newTable = PRM.IMPROPER(data=PRM.readSection(filename, "IMPROPER"))  
            super().__init__(data=self.append(newTable, ignore_index=True).astype({
                     'Kpsi' :float,
                     'psi0' : float
                    }))


#==============================================================================
class NAMDdataEsception(Exception):
    pass
        

#==============================================================================
class NAMDdata:
    ''' Groups PDB, PRM and PSF objects.'''

    def __init__(self, *files):
        self.pdb = PDB()
        self.psf = PSF()
        self.prm = PRM()
        self.pbc = PBC()
        self.network = None
        
        if files:
            self.readFiles(*files)
    
    def readFiles(self, *files):
       
        if len(files) == 0:
            raise NAMDdataEsception("no files given to readFiles function")
        elif len(files) > 4:
            raise NAMDdataEsception("too many files given to readFiles function")
        else:
            for f in files:
                if   ".pdb" in f:
                    self.pdb.readFile(f)
                elif ".psf" in f:
                    self.psf.readFile(f)
                elif ".prm" in f:
                    self.prm.readFile(f)
                elif ".xsc" in f:
                    self.pbc.readFile(f)
                else: 
                    print("file:" + f + "does not have pdb, psf or prm as an extension")


    def loadWolffia(self, wolffia):
        '''
        Converts a Wolffia mixture to a NAMDData self object.

        @param: mix a Wolffia Mixture object.
        '''
        import tempfile

        #print("NAMDdata.loadWolffia() writing files...")
        tdir = tempfile.TemporaryDirectory(suffix="granules_")
        filesSufix = tdir.name + "/namdTempFile"
        wolffia.writeFiles(filesSufix)
        wolffia.writeFiles("namdTempFile")

        if wolffia.hasBox():
            self.readFiles(filesSufix + ".pdb", filesSufix + ".psf", filesSufix + ".prm", filesSufix + ".xsc")
        else:
            self.readFiles(filesSufix + ".pdb", filesSufix + ".psf", filesSufix + ".prm")

        tdir.cleanup()

        return self

 
    def select(self, selection, reverse=False):
        ''' selects atoms acording to selection
        selection is a dictionary with nams of columns of class PDB as keys and 
        lists of values to be selected as values.
        '''

        self.pdb.select(selection, reverse=reverse)
        a_ids = self.pdb.ID
        self.psf.select(a_ids)
        #self.prm.select(Name=Name)
        #self.pbc.select(Name=Name)
            
    
    def get_cycles_oflength(self,n):
        def DFS(graph, marked, n, vert, start, count, cycle): 
            global cycles_sets
            
            marked[vert] = True  # mark the vertex vert as visited 
            cycle.append(vert)
            
            if n == 0:  # if the path of length (n-1) is found 
                if graph.has_edge(vert,start): # Check if vertex vert ends with vertex start 
                    count = count + 1
                    cycles_sets.add(frozenset(cycle))
                    #print(cycle)
        
                # mark vert as un-visited to make it usable again. 
                marked[vert] = False
                cycle.pop()
        
                return count 
            # For searching every possible path of length (n-1) 
            for i in graph.nodes(): 
                if not marked[i] and graph.has_edge(vert,i): 
                    # DFS for searching path by decreasing length by 1 
                    count = DFS(graph, marked, n-1, i, start, count,cycle) 
        
            # marking vert as unvisited to make it 
            # usable again. 
            marked[vert] = False
            cycle.pop()
            return count
        # Counts cycles of length 
        # N in an undirected 
        # and connected graph. 
        def countCycles( graph, n): 
            # all vertex are marked un-visited initially. 
            marked = [False] * (len(graph) + 1)
        
            # Searching for cycle by using v-n+1 vertices s
            count = 0
            for atom in graph.nodes():
                #print("start: " , atom, " con adyacencias ", graph.neighbors(n))
                count = DFS(graph, marked, n-1, atom, atom, count, []) 
                
                # ith vertex is marked as visited and 
                # will not be visited again. 
                #marked[atom] = True
            
            return int(count/2) 
        cycles_sets = set()
        pairs = zip(self.psf.bonds['atom1'], self.psf.bonds['atom2'])
        g = nx.Graph(pairs)
        ngon_bonds = countCycles(g,n)
        polygons = [Polygon(*(self.pdb[self.pdb['ID'].isin(verts)]['ID',['x','y','z']])) for verts in cycle_sets]
        return polygons
        
    def n_gon_connections(self, n):
        pgons = get_cycles_oflength(n)
        return Graph([(p1,p2) for p1 in pgons for p2 in pgons if p1.isneighbor(p2)])


#=============================================================================
if __name__ == "__main__":  # tests
    #pdb = PDB()
    #pdb.readFile("5brr_autopsf.pdb")
    #print(pdb.dtypes)

    #psf = PSF()
    #psf.readFile("5brr_autopsf.psf")
    #print(psf.bonds)

    prm = PRM()
    prm.readFile("par_all27_prot_lipid.prm")
    print(prm.dihedrals)

    #ch = NAMDdata()
    #ch.readFiles("2rvd_autopsf.pdb", "2rvd_autopsf.psf", "par_all36_prot.prm")


