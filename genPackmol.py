#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 11:43:53 2019

@author: nvthaomy
"""
import sys
import numpy
def packmol_noSolvent(f,N,np,pdbList):
    mixturePdb = []
    inpFile = []
    #convert to angstrom
    x,y,z = ((x-.2)*10,(y-.2)*10,(z-.2)*10)
    np = [int(i) for i in np]
    for i,charge in enumerate(f):
        for k,num_poly in enumerate(np):    
            inp = open('mixAH{}_f{}_np{}.inp'.format(N,charge,num_poly),'w')
            inp.write('tolerance 2.0')
            inp.write('\nfiletype pdb')
            inp.write('\noutput AH{}_f{}_np{}.pdb'.format(N,charge,num_poly))
            mixturePdb.append('AH{}_f{}_np{}.pdb'.format(N,charge,num_poly))
            inp.write('\nadd_amber_ter')
            inp.write('\n')
            inp.write('\nstructure {}'.format(pdbList[i]))
            inp.write('\n\tnumber {}'.format(int(num_poly)))
            inp.write('\n\tresnumbers 2')
            inp.write('\n\tinside box 5 5 5 450 450 450')
            inp.write('\nend structure')
            inpFile.append('mixAH{}_f{}_np{}.inp'.format(N,charge,num_poly))
    #writing job file to pack molecules
    with open('packmol.job','w') as job:
        job.write('#!/bin/bash')
        job.write('\n#')
        job.write('\n#$ -V')
        job.write('\n#$ -cwd')
        job.write('\n#$ -j y')
        job.write('\n#$ -S /bin/bash')
        job.write('\n#$ -N packmol')
        job.write('\nmodule load packmol')
        for inp in inpFile:
            job.write('\npackmol < {}'.format(inp))      
    return mixturePdb
def packmol(f,N,np,nw,watermodel,pdbList,x,y,z,ns,s, saltPdbDicts = {'Na+':'na.pdb', 'Cl-':'cl.pdb'}):
    #Calculating weight fraction from nw and np    
    w = ((N*57*np)/(N*57*np+18*nw))
    w = [round(i,2) for i in w]
    inpFile=[]
    mixturePdb = []
    x,y,z = ((x-.2)*10,(y-.2)*10,(z-.2)*10) 
    for i,charge in enumerate(f):
        for k,wtfrac in enumerate(w):    
            if ns[k] > 0.:
                inpName = 'mixAH{}_f{}_{}_w{}_{}{}{}.inp'.format(N,charge,watermodel,wtfrac,int(ns[k]),s[0],s[1])
                pdbName = 'AH{}_f{}_{}_w{}_{}{}{}.pdb'.format(N,charge,watermodel,wtfrac,int(ns[k]),s[0],s[1])
            else:
                inpName = 'mixAH{}_f{}_{}_w{}.inp'.format(N,charge,watermodel,wtfrac)
                pdbName = 'AH{}_f{}_{}_w{}.pdb'.format(N,charge,watermodel,wtfrac)
            mixturePdb.append(pdbName)
            inp = open(inpName,'w')
            inp.write('tolerance 2.0')
            inp.write('\nfiletype pdb')
            inp.write('\noutput {}'.format(pdbName))
            inp.write('\nadd_amber_ter')
            inp.write('\n')
            inp.write('\nstructure {}'.format(pdbList[i]))
            inp.write('\n\tnumber {}'.format(int(np[k])))
            inp.write('\n\tresnumbers 2')
            inp.write('\n\tinside box 3 3 3 {} {} {}'.format(x,y,z))
            inp.write('\nend structure')
            inp.write('\n')
            inp.write('\nstructure {}.pdb'.format(watermodel))
            inp.write('\n\tnumber {}'.format(int(nw[k])))
            inp.write('\n\tresnumbers 2')
            inp.write('\n\tinside box 2 2 2 {} {} {}'.format(x+1,y+1,z+1))
            inp.write('\nend structure')
            if ns[k] > 0.:
                lines = ""
                for ion in s:
                    lines += "\nstructure {ionPdb}".format(ionPdb = saltPdbDicts[ion])
                    lines += "\n\tnumber {n_ion}".format(n_ion = int(ns[k]))
                    lines += "\n\tresnumbers 2\n\tinside box 3 3 3 {x} {y} {z}\nend structure\n".format( x=x, y=y, z=z)
                inp.write(lines)
            inpFile.append(inpName)
    #writing job file to pack molecules
    with open('packmol.job','w') as job:
        job.write('#!/bin/bash')
        job.write('\n#')
        job.write('\n#$ -V')
        job.write('\n#$ -cwd')
        job.write('\n#$ -j y')
        job.write('\n#$ -S /bin/bash')
        job.write('\n#$ -N packmol')
        job.write('\nmodule load packmol')
        for inp in inpFile:
            job.write('\npackmol < {}'.format(inp))      
    return mixturePdb,w        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",nargs='+', type=float, required=True,
                        help="list of ionization fractions e.g. 0 0.2 0.5 1")
    parser.add_argument("-N",type=int, required=True,
                        help="degree of polymerization")
    parser.add_argument("-np",nargs='+', type=float, 
                        help="list of number of PAH chains, must have same length as nw")
    parser.add_argument("-nw", nargs='+', type=float, 
                        help="list of number of water molecules")
    parser.add_argument("-w","--watermodel",
                        help="Water model for simulation (opc,tip3p,spce), default = opc")
    parser.add_argument("-pdb",nargs='+',required=True,
                        help="list of single chain pdb files")
    
    args = parser.parse_args()
    if args.watermodel:
        watermodel = args.watermodel
    else:
        watermodel = 'opc'
    f = args.f
    N = args.N    
    np= numpy.array(args.np)
    nw= numpy.array(args.nw)
    pdbList = args.pdb
    if not len(np) == len(nw) and len(nw) != 0:
        sys.stderr.write('\n np and nw do not have the same number of arguments!')
    newCharge = []
    for index,charge_frac in enumerate(f):
            if charge_frac == 0 or charge_frac == 1:
                charge_frac = int(charge_frac)
            else:
                charge_frac = round(charge_frac,2)
            newCharge.append(charge_frac)
    f = newCharge
    if nw:
        packmol(f,N,np,nw,watermodel,pdbList)
    else:
        packmol_noSolvent(f,N,np,pdbList)
    
    
    
