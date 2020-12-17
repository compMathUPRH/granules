#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 22:09:24 2020

@author: lemuel
"""

#dcdReader for granules

from struct import *
from granules.structure.LAMMPSdata import *

INT_SIZE    = calcsize('i')
FLOAT_SIZE  = calcsize('f')
DOUBLE_SIZE = calcsize('d')


class DCDReader(object):#pq tiene object?

	def __init__(self, filename):
		self.open(filename)
		self.filename = filename
		
	def close(self):
		self.dcdfile.close()
		self.dcdfile = None
	
	def isOpen(self): 
		try: return self.dcdfile != None
		except: return False
		
	def open(self, filename=None):
		if filename == None: filename = self.filename
		if self.isOpen(): self.close()
		self.dcdfile = open(filename, 'rb')
		self.readHeaders()
		self.currentFrame = 0
		
	def restart(self):
		if self.isOpen(): self.close()
		self.open(self.filename)

	def readHeaders(self):
		# == HEADER block ==
		assert(unpack('i', self.dcdfile.read(INT_SIZE))[0] == 84    ) # first integer is 84
		#print(self.dcdfile.read(4))
		assert(self.dcdfile.read(4)                        == b'CORD') #next 4 bytes are 'CORD'

		self.nset    = unpack('i', self.dcdfile.read(INT_SIZE))[0]    # num of sets of coords
		self.istart  = unpack('i', self.dcdfile.read(INT_SIZE))[0]    # starting timestep
		self.stepInt = unpack('i', self.dcdfile.read(INT_SIZE))[0]    # num states between dcd saves
		self.dcdfile . read  (5*INT_SIZE)                        # skip 5 integers
		self.namnf   = unpack('i', self.dcdfile.read(INT_SIZE))[0]    # num of non free atoms
		self.ts      = unpack('f', self.dcdfile.read(FLOAT_SIZE))[0]  # timestep (ps)
		self.wPBC    = unpack('i', self.dcdfile.read(INT_SIZE))[0] == 1  # TRUE if PBCs
		self.dcdfile . read  (9*INT_SIZE)                        # skip 9 integers

		assert(unpack('i', self.dcdfile.read(INT_SIZE))[0] == 84)# should be 84
		#print(self.nset,self.istart,self.stepInt,self.namnf,self.ts,self.wPBC)

		# == TITLE block ==
		tsize = unpack('i', self.dcdfile.read(INT_SIZE))[0]      # size of title block
		assert((tsize-4) % 80 == 0)                         # must be a multiple of 80

		tlines = unpack('i', self.dcdfile.read(INT_SIZE))[0]	    # NUM LINES IN TITLE
		self.title = [self.dcdfile.read(80) for i in range(tlines)]
		self.dcdfile . read  (INT_SIZE)                          # ending size

		#print(self.title)

		assert(unpack('i', self.dcdfile.read(INT_SIZE))[0] == 4) # should be 4
		self.numAtoms = unpack('i', self.dcdfile.read(INT_SIZE))[0]   # num atoms
		self.dcdfile . read  (INT_SIZE)                          # size of arrays

		assert self.namnf == 0  
		#if self.namnf != 0:   # Not ready yet nor needed.  Always namf = 0 in NAND
		#	self.fixAtomInd = [unpack('i', self.dcdfile.read(INT_SIZE))[0] for i in range(self.numAtoms-self.namnf)]
		#	print "fixAtomInd=", self.fixAtomIndDCDReader


	def __iter__(self): return self


	def __next__(self):
		if self.currentFrame == self.nset:
			self.dcdfile.close()
			raise StopIteration

		self.currentFrame += 1

		if self.wPBC:
			self.dcdfile . read  (INT_SIZE)                          # size of arrays
			pbc = [unpack('d', self.dcdfile.read(DOUBLE_SIZE))[0] for i in range(6)]
			self.dcdfile . read  (INT_SIZE)                          # size of arrays
			print("PBC=",pbc)
		else: pbc = None

		if self.currentFrame == 1 or self.namnf == 0:  # always True in NAMD
			freeAtoms = self.numAtoms

			self.dcdfile . read  (INT_SIZE)                          # size of arrays
			X = [unpack('f', self.dcdfile.read(FLOAT_SIZE))[0] for i in range(freeAtoms)] 
			self.dcdfile . read  (INT_SIZE)                          # size of arrays
			#print X

			self.dcdfile . read  (INT_SIZE)                          # size of arrays
			Y = [unpack('f', self.dcdfile.read(FLOAT_SIZE))[0] for i in range(freeAtoms)] 
			self.dcdfile . read  (INT_SIZE)                          # size of arrays
			#print Y

			self.dcdfile . read  (INT_SIZE)                          # size of arrays
			Z = [unpack('f', self.dcdfile.read(FLOAT_SIZE))[0] for i in range(freeAtoms)] 
			self.dcdfile . read  (INT_SIZE)                          # size of arrays
			#print Z

			#print unpack('f', dcdfile.read(FLOAT_SIZE))[0]
			return (X,Y,Z,pbc)
 
class DCDTrajectoryReader(DCDReader):
	#from wolffia.interface.main.History import NanoCADState
	def __init__(self, simulation, filename):
		super(DCDTrajectoryReader, self).__init__(filename)
		self.simulation = simulation

class TrayectoryReader(DCDTrajectoryReader):
	pass # eventually will work for other formats than DCD

class Trajectory(DCDTrajectoryReader):
	def __init__(self, filename):
		super(DCDTrajectoryReader, self).__init__(filename)
		self.filename = filename
	'''
	def __getitem__(self,a, b):
		dequeue(None,)
		for X,Y,Z in r:
	'''
			
	def head(self, n=1):
		''' Returns a list of the first n frames of the trajectory.
		'''
		i=0
		result=[]
		for X,Y,Z,pcb in r:
			result.append(zip(X,Y,Z,pcb))
			i += 1
			if i == n: break
		return result
		
	def tail(self, n=1):
		''' Returns a list of the first n frames of the trajectory.
		'''
		from collections import deque
		self.restart()  # eventually this will be handled by "seek"  operation
		result=deque([],n)
		for X,Y,Z in r:
			result.append(zip(X,Y,Z))
		return list(result)            
    
if __name__ == '__main__':
    ################TEST######################
    print("Corrida#1")
    r = Trajectory("chignolin.dcd")
    #r = DCDReader("chignolin.dcd")
    l = LammpsData()
    '''
    for frame in r:
        l.atomproperty.atoms['x'] = frame[0]
        l.atomproperty.atoms['y'] = frame[1]
        l.atomproperty.atoms['z'] = frame[2]
        
    '''    
    print("Corrio")
            
    
    print("=========TODOS================")
    i = 0
    #for X,Y,Z,pbc in r:
    for tuples in r.head(10):
        #print X,Y,Z
        #tuples = zip(X,Y,Z)
        print(tuples)
        i += 1
    		#print [item for sublist in tuples for item in tuples]
    print(i)
    '''
    print("==============================")
    h=r.tail(3)
    for c in h:
        print(c[:2])
           
    '''        
            
        
        