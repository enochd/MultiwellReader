#/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script contains a number of functions useful in generating plots and fits from 
MW output
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from scipy.constants import physical_constants

#def loadData(species, path):
def loadData(path, species):
"""
From a MW output file of a pdep calculation located at 'path'
plot the 'species' concn as a fucntion of time and extrapolate 
the phenomenological rate coefficient. 
"""
    try:
        f = open(path, 'r')
    except IndexError:
        print('File not found')
    
    f = open(path, 'r')
    line = f.readline()
    while line != '':
        # The data we want is in the end of the output
        if '                Name:' in line:
            #i = 0
            for i in range(len(line.split())):
                if species == line.split()[i]:
                    concline = i 
                    #i += 1
            
        if 'TEMPERATURES' in line:
            line = f.readline() 
            T = float(line.split()[3])
        if 'Press' in line: 
            P = float(line.split()[2])
        if 'CPU :' in line: 
            cpu = float(line.split()[5])       
        if 'Collisns' in line:     #purposefully a typo, fyi!
            print 'found a sucessfully completed MW Job'
            time = []
            conc = []
            line = f.readline()
            line = f.readline()
            break
        line = f.readline() 
    while line != '':
        if float(line.split()[concline*3-1]) == 0.0:
            print 'uh oh'
            break
        time.append(float(line.split()[0]))
        conc.append(float(line.split()[concline*3-1]))
        # Read the next line in the file  
        line = f.readline()    
    f.close()
    
    Conc = np.array(conc)
    Time = np.array(time) 
    #Test the fitting:
    #time = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    #conc = [110132.328974034, 742.0657955129, 140.1581244726, 60.9124698035, 36.9452804947, 26.4724502524, 20.863669418, 17.4517147873, 15.1886588876, 13.5914091423, 12.4103254231, 11.5048794545, 10.7905276697, 10.2136353513, 9.7386702053, 9.3412297872, 9.0040385688, 8.7145449932, 8.4634230016, 8.2436063535]
    #Time = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
    #Conc = np.array([110132.328974034, 742.0657955129, 140.1581244726, 60.9124698035, 36.9452804947, 26.4724502524, 20.863669418, 17.4517147873, 15.1886588876, 13.5914091423, 12.4103254231, 11.5048794545, 10.7905276697, 10.2136353513, 9.7386702053, 9.3412297872, 9.0040385688, 8.7145449932, 8.4634230016, 8.2436063535])

    weights=None
    threeParams=True
    T0=1
    # kunits = string
    if threeParams:
        A = np.zeros((len(Time),2), np.float64)
        A[:,0] = np.ones_like(Time)
        A[:,1] = Time                 # natural log
        b = np.log(Conc)    #conc cant go to zero! 
        x, residues, rank, s = np.linalg.lstsq(A,b)
        '''
        x = least-squares solution
        b - array of dependent variables
        s - singular values of A
        '''
        r2 = 1 - residues / (b.size * b.var())
        intercept = np.exp(x[0]) 
        k = x[1]#/np.log(Tlist)
        print x
        print "A, k,are " + str(intercept) + " " + str(k) 
        print "R^2 is " + str(r2)  
        print "s is " + str(s) 
        print "MW took " + str(cpu) + 'hrs'
        y = intercept*np.exp(k * Time)

    #TO DO LIST:   
    # remove first few points (where steady state conditions have not yet been achieved)
    # save branching to output  
    fig = plt.figure()
    fig.suptitle(str(species) + ' decomposition', fontsize=20, fontweight='bold', style='italic')
    ax = fig.add_subplot(111)    
    fig.subplots_adjust(top=0.85)
    #this will fail if arrays are not large enough:
    if len(time) > 70 and len(conc) > 50:
        ax.text(time[70], conc[50], r'$R^2$ = '+"%1.4f" % float(r2), fontsize=15)
    else:
        ax.text(time[1], conc[1], r'$R^2$ = '+"%1.4f" % float(r2), fontsize=15)         
    plt.semilogy(Time, y, 'r', linewidth=2.0)
    plt.semilogy(time, conc, 'yo-')
    plt.title('T = ' + "%5.0f" % float(T) + 'K; P = ' + str(P) + ' atm;' + ' k is ' + "%1.4e" % -k + ' ' + r'$s^-$' + r'$^1$')
    plt.ylabel('concentration', fontsize=20)
    plt.xlabel('time, s', fontsize=20)
    #ax.text(0, 0.2, str(-k), fontsize=10)
    ax.text(time[2], conc[1]/2, 'CPU time: ' + str(cpu) + ' hrs', fontsize=15)
    fname = str(species) + '_' + str(P) + '_' + str(T) + '.pdf'
    plt.savefig(fname,format='pdf') 
    plt.show() #must come AFTER savefig
    
    print fname
    return T, P, k            


def loadAllData(path, Temp, species):
    """
    read an input file for the path to MW out files. 
    read species of interest
    extract k(T,P), if possible, and save rate info to output
    """   

    try:
        f = open(path, 'r')
    except IndexError:
        print('File not found')

    f = open(path, 'r')
    print path 
    line = f.readline() 
    while line != '':
        print 'blah'
        # open and process each output file
        pres = line.split()[0]
        #pres = float(line.split()[0])
        print 'proccessing '+ str(pres) + ' atm'
        path2 = path[:-14] + str(Temp) + '/' + str(pres) + '.out'
        print path2
        loadData(path2, species) 
        
        line = f.readline()
         
    f.close()    
    return 

#loadData('/home/enoch/MW_comparisons_cantherm/ethoxy/hydroxyethyl/RRHO/900/0.01.out', 'ch3choh')
#loadData('/home/enoch/MW_comparisons_cantherm/ethoxy/RRHO/800/0.001.out', 'ch3ch2o') 
#loadData('/home/enoch/MW_comparisons_cantherm/Methoxy_CH2OH/methoxy_single_channel/pharos_imports/678K_RRHO/100.out', 'ch3o') 
loadData('/home/enoch/multiwell2014/CPD+CPDyl/coll.out', 'out')
loadData('/home/enoch/multiwell2014/CPD+CPDyl/coll.out', 'w3')
#loadAllData('/home/enoch/MW_comparisons_cantherm/ethoxy/hydroxyethyl/RRHO/inp_pylist.txt', 900, 'ch3choh')
#loadAllData('/home/enoch/MW_comparisons_cantherm/ethoxy/RRHO/inp_pylist.txt', 1000, 'ch3ch2o')

def fitToMwData(conc, time):
    Conc = np.array(conc)
    Time = np.array(time)    
    
    pass

def obtainBranching(conc, time):
    Conc = np.array(conc)
    Time = np.array(time)    
    
    pass

def getRateCoefficient(Tlist,Plist,kexplist):
    """
    Test the MultiPDepArrhenius.getRateCoefficient() method.
    """
    for i in range(Tlist.shape[0]):
        for j in range(Plist.shape[0]):
            kexp = kexplist[i,j]
#            return kexp

pass
'''
Tlist = np.array([200,400,600,800,1000,1200,1400,1600,1800,2000]) # b = np.arange(1, 9, 2), start, end (exclusive), step
Plist = np.array([1e4,1e5,1e6])
kexplist = np.array([
[2.85400e-08, 4.00384e-03, 2.73563e-01, 8.50699e+00, 1.20181e+02, 7.56312e+02, 2.84724e+03, 7.71702e+03, 1.67743e+04, 3.12290e+04],
[2.85400e-07, 4.00384e-02, 2.73563e+00, 8.50699e+01, 1.20181e+03, 7.56312e+03, 2.84724e+04, 7.71702e+04, 1.67743e+05, 3.12290e+05],
[2.85400e-06, 4.00384e-01, 2.73563e+01, 8.50699e+02, 1.20181e+04, 7.56312e+04, 2.84724e+05, 7.71702e+05, 1.67743e+06, 3.12290e+06],
]).T
getRateCoefficient(Tlist,Plist,kexplist)
'''
#to do: make python read in Tlist and klist
Tlist = np.array([800,900,1000,1200,1500])
klist = np.array([1.66E+05, 1.10E+06, 5.10E+06, 5.27E+07, 5.68E+08])
#klist = np.array([2.85400e-08, 4.00384e-03, 2.73563e-01, 8.50699e+00, 1.20181e+02, 7.56312e+02, 2.84724e+03, 7.71702e+03, 1.67743e+04, 3.12290e+04])
#fitToData(Tlist,klist)
#print np.size(klist)
#print np.size(Tlist)




def fitToData(Tlist, klist):
    """
    Fit the Arrhenius parameters to a set of rate coefficient data `klist`
    in units of `kunits` corresponding to a set of temperatures `Tlist` in 
    K. A linear least-squares fit is used, which guarantees that the 
    resulting parameters provide the best possible approximation to the 
    data.
    """
    weights=None
    threeParams=True
    T0=1
#        kunits = string
    if threeParams:
        A = np.zeros((len(Tlist),3), np.float64)
        A[:,0] = np.ones_like(Tlist)
        A[:,1] = np.log(Tlist / T0)                 # natural log
        A[:,2] = -1.0 / 1.9872 / Tlist
        b = np.log(klist)
#        A[:,2] = -1 / scipy.constants.physical_constants(R) / Tlist
#        A[:,2] = -1.0 / scipy.constants.R / Tlist
        # Perform a least squares fit to the user given rate data
        x, residues, rank, s = np.linalg.lstsq(A,b)
        '''
        x = least-squares solution
        b - array of dependent variables
        s - singular values of A
        '''
        Apre = np.exp(x[0]) 
        n=x[1]#/np.log(Tlist)
        Ea=x[2]
        print "A, n, and Ea are " + str(Apre) + " " + str(n) + " " + str(Ea)
        # Determine covarianace matrix to obtain parameter uncertainties
        count = klist.size
        cov = residues[0] / (count - 3) * np.linalg.inv(np.dot(A.T, A))
        t = scipy.stats.t.ppf(0.975, count - 3)
        print A
        plt.semilogy(1000./Tlist, klist,'o', label = 'original data)')
        y= Apre*(Tlist**n)*np.exp(-Ea/1.9872/Tlist)
        print y
        plt.semilogy(1000./Tlist,y)
        y2 = 7E-8*Tlist**6.3*np.exp(-19750/1.9872/Tlist)
        print y2
        plt.semilogy(1000./Tlist,y2)
        plt.grid(True)
    return    b, x, Apre, n , Ea, y, y2

