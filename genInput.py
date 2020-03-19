#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:18:18 2019

@author: nvthaomy
"""
import numpy
import sys
import os
import writeOpenMMinput_amber,writeOpenMMinput_cgen, buildPAH
import time

if __name__ == "__main__":
    '''Need to be executed in a folder with 
       PAH.lib, pdb of water, pdbs of up and down AH monomer (both ionized and neutral verions) if a single chain need to be built
       charm36.ff and ions.mdp if using cgen'''
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",nargs='+', required=True,
                        help="list of deprotonation fractions e.g. 0 0.2 0.5 1, if pdb of\
                        a single chain already exists, need to enter the same decimal place\
                        as the file name")
    parser.add_argument("-N",type=int, required=True,
                        help="degree of polymerization")
    parser.add_argument("-Pm","--Pm",type=float, default = 1, 
                        help="Meso diad fraction, isotactic if = 1, syndiotactic if = 0")
    parser.add_argument("-r","--random", action = "store_true",
                        help="Random deprotonation, default pattern is random")
    parser.add_argument("-e","--evendist", action = "store_true",
                        help="Evenly distributed deprotonation")
    parser.add_argument("-w","--watermodel", default = 'opc',
                        help="Water model for simulation (opc,tip3p,spce)")
    parser.add_argument('-t','--temperature',type=float, nargs='+', required = True,
                        help = "list of temperature in K at which simulation is run")
    parser.add_argument("-np",nargs='+', type=float, required=True,
                        help="list of number of PAH chains, must have same length as nw")
    parser.add_argument("-nw", nargs='+', type=float, required=True,
                        help="list of number of water molecules")
    parser.add_argument("-ff", default = "ff99", choices = ['gaff2','ff99'],
                        help="force field")
    parser.add_argument("-l", nargs='+',
                        help="can have multiple libraries, tleap library for PAH monomers: PAH, PAH_avg,PAH1,etc.")
    parser.add_argument("-xyz",nargs ="+",type = float,
                        help="periodic box size in nm")
    parser.add_argument("-a", action = "store_true",
                        help="if use anisotropic MC barostat, acts on the z direction")
    parser.add_argument("-ns", nargs='+', type=float, default = [0],
                        help="list of number of salt molecules (cation and anion pair)")
    parser.add_argument("-s", nargs=2, type=str, default = ['Na+', 'Cl-'],
                        help="list of cation and anion")

    args = parser.parse_args() 
    
    ff = args.ff
    watermodel = args.watermodel
    np=args.np
    np= numpy.array(args.np)
    nw= numpy.array(args.nw)
    f = args.f
    N = args.N
    Pm = args.Pm
    T = args.temperature
    PAHLib = args.l
    ns = args.ns
    s = args.s
    singleChainPdb = []
    pattern = 'random' #random deprotonation
    if args.a:
        anisoP = True
    else:
        anisoP = False

    if args.evendist:
        pattern = 'even'
    if not len(np) == len(nw):
        sys.stderr.write('\n np and nw do not have the same number of arguments!')
    #weight fraction from MW of neutral monomers

    #estimate box size from water density if not provided:
    if not args.xyz:
        ntot = nw[0] + N*np[0]
        vol = float(ntot)/33.36
        x = vol**(1./3.)
        y = x
        z = x
    else:
        x,y,z = (args.xyz[0],args.xyz[1],args.xyz[2])

     
    w = ((N*57*np)/(N*57*np+18*nw)) #calculate weight fraction
    w = [round(i,2) for i in w]
    print ('wt fraction {}'.format(w)) 
    needToBuildSingleChain = raw_input('\nNeed to make single chain pdb?  (y/n) ')
    if needToBuildSingleChain == 'y':
        print "\nWriting build.in file to construct polymer chains in tleap"
        f = [float(i) for i in f]
        build_file,singleChainPdb,f = buildPAH.buildPAH(f,N,Pm,pattern)
        while not os.path.exists(build_file):
            time.sleep(1)
        print('\nBuilding polymer in tleap ...')
        os.system('\ntleap -s -f {} > build_AH{}.out'.format(build_file,N)) #call tleap to build PAH
    else:
        singleChainPdb = raw_input('\nEnter list of single chain pdbs (e.g.:AH22_f0.pdb AH22_f1.pdb): ')
        singleChainPdb = list(singleChainPdb.split())
        for pdb in singleChainPdb:
            if not os.path.exists(pdb):
                print "\n{} does not exist!".format(pdb)
                break
        #for charge_frac in f:
        #    singleChainPdb.append('AH{}_f{}.pdb'.format(N,charge_frac))
        f = [float(i) for i in f]
        newCharge = []
        for charge_frac in f:
            if abs(charge_frac - 0) < 10**(-2) or abs(charge_frac - 1) < 10**(-2): 
                charge_frac = int(charge_frac)
            else:
                charge_frac = round(charge_frac,2)
            newCharge.append(charge_frac)
        f = newCharge
    if ff in ['gaff2','ff99']:
        writeOpenMMinput_amber.main(f,N,np,nw,watermodel,singleChainPdb,T,PAHLib,ff, x,y,z, anisoP, ns, s)
    elif ff == 'cgen':
        writeOpenMMinput_cgen.main(f,N,np,nw,w,watermodel,singleChainPdb,T)
        
