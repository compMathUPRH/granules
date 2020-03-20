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
    ''' Pandas DataFrame that stores data in PDB format (only the ATOMS section).
    '''

    def __init__(self,data=None, dtype=None, copy=False):
        super(PDB, self).__init__(data=data, 
                                  copy=copy,
                                  columns=['RecName','ID','Name','AltLoc','ResName',
                                           'ChainID','ResSeq','iCode','x','y','z',
                                           'Occupancy','TempFactor','Element','Charge']
        )

    def readFile(self, filename):
        ''' reads PDB file and appends to self

            Parameter
            ----------
            filename : str
                name of file
        '''

        #Abrir el archivo para leer datos
        arch = open(filename, 'r')
        data = []
        
        ranges = [6,6,4,1,3,2,4,4,8,8,8,6,6,11,2]
        
        for linea in arch:
            
            #Solo queremos la informacion
            if linea[:4] == 'ATOM' or linea[:4] == 'HETA': 
                
                #Recoger la informacion mediante su posicion
                info = []
                start = 0
                stop = 0
                        

                for i in ranges:
                    stop += i
                    #print('\nEn posicion',start,':',stop,'se encuentra ->>')
                    l = linea[start:stop].strip()
                    if l == '' or l == '<0>': info.append(np.nan)
                    else: info.append(l)
                    start += i

                #Cada vez anadir nueva informacion a los datos
                #print(info)
                data.append(info)
                       
        arch.close() 
        newTable = PDB(data=data)  

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




class PSF:
    ''' Pandas DataFrames that store data in PSF format.
        See specification in http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html
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
        self.atoms.readSection(filename)
        self.bonds.readSection(filename)
        self.angles.readSection(filename)
        self.dihedrals.readSection(filename)
        self.impropers.readSection(filename)
        self.cross_terms.readSection(filename)

    @staticmethod
    def readSection(filename, section, tupleLength, itemsPerLine):
        ''' reads section 'section' PSF file

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
        flag = True
        

        # find desired section
        for linea in arch:
            if '!N'+section in linea: break
        #print("\n\n\nlinea encontrada: ", linea)


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


    class ATOM(pd.DataFrame):
        ''' ATOM section of the PSF '''

        def __init__(self,data=None, dtype=None, copy=False):
            super(PSF.ATOM, self).__init__(data=data, copy=copy, columns=[
                        'ID','RecName','ChainID', 'ResName', 'Name', 
                        'Type', 'Charge', 'Mass', 'Unused'])

        def readSection(self, filename):
            newTable = PSF.ATOM(data=PSF.readSection(filename, "ATOM", 9, 1))  
            #print(newTable)
            # set column types and append to existing table
            super().__init__(data=self.append(newTable, ignore_index=True).astype({
                     'ID'     :int,
                     'Charge' :float,
                     'Mass'   :float
                    }))


    class BOND(pd.DataFrame):
        ''' BOND section of the PSF'''

        def __init__(self,data=None, dtype=None, copy=False):
            super(PSF.BOND, self).__init__(data=data, copy=copy, columns=[
                    'atom1','atom2'])

        def readSection(self, filename):
            data=PSF.readSection(filename, "BOND", 2, 4)
            newTable = pd.DataFrame(data).dropna()
            pair1 = newTable[[0,1]]

            pair1.columns=['atom1','atom2']

            super().__init__(data=self.append(pair1, ignore_index=True).astype({
                     'atom1'     :int,
                     'atom2' :int
                    }))

    class THETA(pd.DataFrame):
        ''' THETA section of the PSF'''

        def __init__(self,data=None, dtype=None, copy=False):
            super(PSF.THETA, self).__init__(data=data, copy=copy, columns=[
                    'atom1','atom2','atom3'])

        def readSection(self, filename):
            data=PSF.readSection(filename, "THETA", 3,3)
            newTable = pd.DataFrame(data).dropna()
            pair1 = newTable[[0,1,2]]
            pair1.columns=['atom1','atom2','atom3']

            super().__init__(data=self.append(pair1, ignore_index=True).astype({
                     'atom1'     :int,
                     'atom2'     :int,
                     'atom3' :int
                    }))

    class PHI(pd.DataFrame):
        ''' PHI section of the PSF'''

        def __init__(self,data=None, dtype=None, copy=False):
            super(PSF.PHI, self).__init__(data=data, copy=copy, columns=[
                    'atom1','atom2','atom3','atom4'])

        def readSection(self, filename):
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


    class IMPHI(pd.DataFrame):
        ''' IMPHI section of the PSF'''

        def __init__(self,data=None, dtype=None, copy=False):
            super(PSF.IMPHI, self).__init__(data=data, copy=copy, columns=[
                    'atom1','atom2','atom3','atom4'])

        def readSection(self, filename):
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


    class CRTERM(pd.DataFrame):
        ''' CRTERM section of the PSF'''

        def __init__(self,data=None, dtype=None, copy=False):
            super(PSF.CRTERM, self).__init__(data=data, copy=copy, columns=[
                    'atom1','atom2','atom3','atom4'])

        def readSection(self, filename):
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


class PRM:
    ''' NAMDdata force field from a PRM file.
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

        self.bonds.readSection(filename)
        self.angles.readSection(filename)
        self.dihedrals.readSection(filename)
        self.impropers.readSection(filename)
        self.nonbonded.readSection(filename)

   
    @staticmethod
    def readSection(filename, section):
        ''' reads section 'section' PRM file

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
        flag = True
       

        # find desired section
        '''
        linea = arch.readline().split("!")[0].strip()
        while not section == linea[:len(section)]: 
           linea = arch.readline().split("!")[0].strip()
        '''
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
        #print(data)
        return data


    class NONBONDED(pd.DataFrame):
        def __init__(self,data=None, dtype=None, copy=False):
            super(PRM.NONBONDED, self).__init__(data=data, copy=copy, columns=[
                'Type','epsilon','Rmin2','epsilon1_4','Rmin2_1_4'])

        def getCoeffs(self):
            # retrieve info from PRM
            prmFF = self.copy()
            prmFF.set_index('Type', inplace=True)
            return prmFF


        def readSection(self, filename):

            data = list()

            #print("Entre a funcion")

            d = PRM.readSection(filename, "NONBONDED")

            for list_ in d:
                #print(list_)
                if len(list_) < 7:
                    data.append([list_[0], -float(list_[2]), list_[3], -float(list_[2]), list_[3]])
                else:
                    data.append([list_[0], -float(list_[2]), list_[3], -float(list_[5]), list_[6]])

            #print(data)
            newTable = PRM.NONBONDED(data=data)  

            # modificar newTable para que todas las listas tengan 4 números
            #print(newTable)

            super().__init__(data=self.append(newTable, ignore_index=True).astype({
                     'epsilon' :float,
                     'Rmin2' :float,
                     'epsilon1_4' :float,
                     'Rmin2_1_4' :float
                    }))




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
            newTable = PRM.BONDS(data=PRM.readSection(filename, "BONDS"))  
            super().__init__(data=self.append(newTable, ignore_index=True).astype({
                     'Kb' :float,
                     'b0'   :float
                    }))


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
            newTable = PRM.DIHEDRALS(data=PRM.readSection(filename, "DIHEDRALS"))  
            super().__init__(data=self.append(newTable, ignore_index=True).astype({
                     'Kchi'  :float,
                     'n'     :int,
                     'delta' : float
                    }))

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
            newTable = PRM.IMPROPER(data=PRM.readSection(filename, "IMPROPER"))  
            super().__init__(data=self.append(newTable, ignore_index=True).astype({
                     'Kpsi' :float,
                     'psi0' : float
                    }))





class NAMDdata:
    ''' Groups PDB, PRM and PSF objects.'''

    def __init__(self, *files):
        self.pdb = PDB()
        self.psf = PSF()
        self.prm = PRM()
        
        if files:
            self.readFiles(*files)
        
    ___commented_old_function = """
    def readFiles(self, file_pdb=None, file_psf=None, file_prm=None):
        #
        if file_pdb != None: self.pdb.readFile(file_pdb)
        if file_psf != None: self.psf.readFile(file_psf)
        if file_prm != None: self.prm.readFile(file_prm)
    """
    
    def readFiles(self, *files):
       
        if len(files) == 0:
            print("no files given to readFiles function")
        elif len(files) > 3:
            print("too many files given to readFiles function")
        else:
            for f in files:
                if ".pdb" in f:
                    self.pdb.readFile(f)
                elif ".psf" in f:
                    self.psf.readFile(f)
                elif ".prm" in f:
                    self.prm.readFile(f)
                else: 
                    print("file:" + f + "does not have pdb, psf or prm as an extension")
        
        


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


