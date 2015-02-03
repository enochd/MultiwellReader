#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
# This script 
# 
#
################################################################################

import math
import numpy
import os.path
from QchemLog import QchemLog
################################################################################


def multiwellDensumMakeFile(path, folderPath):
    '''
    Automatically generate all the DensData files needed for a Multiwell Simulation.

    Given a text file with a list of QChem output file names, search the folder 
    provided by 'path' and make input files for the DenSum utility in Multiwell. 
    Assumes, for now, that the user has performed a Qchem calculation, although this 
    can/will be changed later to incorporate Gaussian jobs. As such, the user must 
    also provide a folderPath where the Qchem output is located (path and folderPath dont
    need to be the same, of course). The Densum input files are generated where this file
    is executed. This program also stays in tune with Multiwell, in that all 
    molecules are treated as symmetric tops. Thus, a 1D external rotor is identified via 
    I_1D,ext != Ib ~ Ic (assuming two external moments of enertia are alike and the third
    is unique)

    The user is prompted to provide some additional information: 
     - vibrational scaling factor for level of theory used in quantum chemistry job, 
    '''
    try:
        f = open(path, 'r')
        print('found the input list text file')
    except IndexError:
        print('Input file not found')    
    try:
        os.path.isdir(folderPath) == True
        print('Found the folder containing the Qchem output files')
    except IndexError:
        print('Qchem output file folder not found')    
    # Prompt user for the scaling factor for the frequencies
    scalingFactor = raw_input("Please enter the scaling factor for the level of theory used: ")    
    scalingFactor = float(scalingFactor)

    for line in f: 
        print 'Making {0} denssum input'.format(line[:-1])
        filePath = folderPath + str(line[:-1]) + '.out'
        fileName = str(line[:-1]) + '.out'
        try:
            f = open(filePath, 'r')
            print 'found the Qchem output files for {0}'.format(fileName)
        except IndexError:
            print 'Qchem output file {0} not found. Are you sure the label in inp.list is correct?'.format(fileName)
        #print filePath
        freq = QchemLog(filePath).loadFreq()
#        print freq
        w = open(line, 'w')
        w.write(line)
        w.write("'"+str(line[:-1])+"'" + '\n')
        #Here I assume that all molecules will be polyatomic nonlinear, for now
        #nDof = QchemLog(filePath).getNumberOfAtoms()*3-6+1  #+1 for the kro
        nDof = len(freq) + 1 # add one degree of freedome for the krotor
        w.write(str(nDof) + ' 1' + " 'HAR'" + " 'CM-1'" + '\n' + '10 300 500 85000'+ '\n')
        for i in range(len(freq)):
            w.write(str(i+1) + ' ' + 'vib ' + str(freq[i]*scalingFactor) + ' 0.0 1' + '\n')
        w.write(str(len(freq)+1) + ' kro ' + "%1.4f" % float(QchemLog(filePath).loadKro()) + ' 0.0 1' + '\n')    
        w.close() 

#path = '/home/enoch/Link to Dropbox (MIT)/vinyl+butadiene/multiwell/inp.list'
#folderPath = '/home/enoch/Link to Dropbox (MIT)/vinyl+butadiene/cantherm/species/'
#multiwellDensumMakeFile(path, folderPath)
multiwellDensumMakeFile('/home/enoch/Link to Dropbox (MIT)/vinyl+butadiene/multiwell/inp.list', '/home/enoch/Link to Dropbox (MIT)/vinyl+butadiene/cantherm/species/')
