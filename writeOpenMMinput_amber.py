#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:12:58 2019

@author: nvthaomy
"""
import sys, glob
import os, shutil
import loadFF
import genPackmol as packmol
import time
def writeJobfile(temp,charge,w,N):
    '''write slurm job file for Pod cluster'''
    name = 'job{}'.format(int(temp))
    charge = ''.join(str(charge).split('.'))
    w = ''.join(str(w).split('.'))
    with open(name,'w') as job:
        job.write('#!/bin/bash')
        job.write('\n#SBATCH -N 1 --partition=gpu --ntasks-per-node=1')   
        job.write('\n#SBATCH --gres=gpu:1')
        job.write('\n#SBATCH -J A{}f{}w{}{}K\n'.format(N,charge,w,int(temp)))
        job.write('\ncd $SLURM_SUBMIT_DIR\n')        
        job.write('\n/bin/hostname')
        job.write('\nsrun --gres=gpu:1 /usr/bin/nvidia-smi')
        job.write('\nexport PATH="/home/mnguyen/anaconda3/bin:$PATH"')
        job.write('\npython sim{}.py'.format(int(temp)))
    return name
def writeOpenMMinput(temp,top,crd,x,y,z,anisoP,np,N):
    name = 'sim{}.py'.format(int(temp))
    with open(name,'w') as sim:
        sim.write('from sys import stdout\nimport numpy as np')
        sim.write('\nimport mdtraj, simtk')
        sim.write('\nfrom simtk.openmm import *')
        sim.write('\nfrom simtk.openmm.app import *')
        sim.write('\nfrom simtk.unit import *\n')
        sim.write('\n# Input Files\n')
        sim.write("\ntopName = '{}'".format(top))
        sim.write("\nTrajFile = 'trajectory{}.dcd'\nThermoLog = 'log{}.txt'".format(int(temp),int(temp)))
        sim.write('\nprmtop = AmberPrmtopFile(\'{}\')'.format(top))
        sim.write('\ninpcrd = AmberInpcrdFile(\'{}\')\n'.format(crd))
        sim.write('\n# System Configuration\n')
        sim.write('\nnonbondedMethod = LJPME')
        sim.write('\nnonbondedCutoff = 1.0*nanometers\ntailCorrection = True')
        sim.write('\newaldErrorTolerance = 0.0001')
        sim.write('\nconstraints = HBonds')
        sim.write('\nrigidWater = True')
        sim.write('\nconstraintTolerance = 0.000001\nbox_vectors = np.diag([{},{},{}]) * nanometer'.format(x,y,z))
        sim.write('\n# Integration Options\n')
        sim.write('\ndt = 0.002*picoseconds')
        sim.write('\ntemperature = {}*kelvin'.format(temp))
        sim.write('\nfriction = 1.0/(100.*dt)')
        sim.write('\npressure = 1.0*atmospheres')
        sim.write('\nbarostatInterval = 25\n')
        sim.write('\n# Simulation Options\n')
        sim.write('\nsteps = 2e7')
        sim.write('\nequilibrationSteps = 100000')
        sim.write('\nplatform = Platform.getPlatformByName(\'CUDA\')')
        sim.write('\nplatformProperties = {\'Precision\': \'mixed\'}')
        sim.write('\ndcdReporter = mdtraj.reporters.DCDReporter(TrajFile, 5000)')
#       sim.write('\ndcdReporter = DCDReporter(\'trajectory{}.dcd\', 5000)'.format(int(temp)))
        sim.write('\ndataReporter = StateDataReporter(ThermoLog, 1000, totalSteps=steps, step=True, speed=True, progress=True, remainingTime=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator=\'\\t\')\n')
        sim.write('\n# Prepare the Simulation\n')
        sim.write('\nprint(\'Building system...\')')
        sim.write('\ntopology = prmtop.topology')
        sim.write('\npositions = inpcrd.positions')
        sim.write('\nsystem = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)')
        if anisoP:
            sim.write('\nsystem.addForce(MonteCarloAnisotropicBarostat(3*[pressure],temperature,False,False,True,barostatInterval))')
        else:
            sim.write('\nsystem.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))')
        sim.write('\nintegrator = LangevinIntegrator(temperature, friction, dt)')
        sim.write('\nintegrator.setConstraintTolerance(constraintTolerance)')
        sim.write('\n#simulation = Simulation(topology, system, integrator, platform)')
        sim.write('\nsimulation = Simulation(topology, system, integrator, platform, platformProperties)')
        sim.write('\nsimulation.context.setPositions(positions)')
        sim.write('\nif inpcrd.boxVectors is not None:')
        sim.write('\n\tsimulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)\n')
        sim.write("""\nsimulation.context.setPeriodicBoxVectors(*box_vectors)
forces = system.getForces()
for force in forces:
    if isinstance(force,simtk.openmm.openmm.NonbondedForce):
        nonbondedforce = force
nonbondedforce.setUseDispersionCorrection(tailCorrection)
nonbondedforce.updateParametersInContext(simulation.context)
forces = system.getForces()
for force in forces:
    if isinstance(force,simtk.openmm.openmm.NonbondedForce):
        nonbondedforce = force
print('getUseDispersionCorrection')
print(nonbondedforce.getUseDispersionCorrection())""")

        sim.write('\n#Restart and Check point')
        sim.write('\nsimulation.saveState(\'output{}.xml\')\n'.format(int(temp)))
        sim.write('\n#to load state: simulation.loadState(\'output.xml\')')
        sim.write('\nsimulation.reporters.append(CheckpointReporter(\'checkpnt{}.chk\', 5000))\n'.format(int(temp)))
        sim.write('\n# Minimize and Equilibrate\n')
        sim.write('\nprint(\'Performing energy minimization...\')')
        sim.write('\nsimulation.minimizeEnergy()')
        sim.write('\nprint(\'Equilibrating...\')')
        sim.write('\nsimulation.context.setVelocitiesToTemperature(temperature)')
        sim.write('\nsimulation.step(equilibrationSteps)\n')
        sim.write('\n# Simulate\n')
        sim.write('\nprint(\'Simulating...\')')
        sim.write('\nsimulation.reporters.append(dcdReporter)')
        sim.write('\nsimulation.reporters.append(dataReporter)')
        sim.write('\nsimulation.currentStep = 0')
        sim.write('\nsimulation.step(steps)')
        sim.write('\nsimulation.saveState(\'output{}.xml\')\n'.format(int(temp)))
        sim.write("""\n
#analyzing
import sys
sys.path.append('/home/mnguyen/bin/PAAMD/')
import analysis

NAtomsPerChain = None #None if analyzing AA
NP = {np}
top = topName
DOP = {N}
fi = 'openmm' #'lammps' or 'openmm' 
analysis.getStats(TrajFile, top, NP, ThermoLog, DOP = DOP, NAtomsPerChain = NAtomsPerChain, StatsFName = 'AllStats.dat',
            RgDatName = 'RgTimeSeries', ReeDatName = 'ReeTimeSeries',RgStatOutName = 'RgReeStats', Ext='.dat',
             fi = fi, obs = [ 'Potential_Energy_(kJ/mole)',   'Kinetic_Energy_(kJ/mole)',   'Temperature_(K)',   'Box_Volume_(nm^3)',   'Density_(g/mL)'], cols = None,
             res0Id = 0, stride = 1, autowarmup = True, warmup = 100, plot = False)
""".format(np=int(np),N=N))
    return name

def main(f,N,np,nw,watermodel,singleChainPdb,T,PAHLib,ff, x,y,z, anisoP, ns, s):
    mixturePdb,wtFraction = packmol.packmol(f,N,np,nw,watermodel,singleChainPdb,x,y,z,ns,s)
    print('\nSubmitting packmol job ...')
    for pdb in singleChainPdb: #checking if single chain pdb files and packmol.job exist
        while not os.path.exists(pdb) or not os.path.exists('packmol.job'):  
            time.sleep(1)
    os.system('qsub packmol.job') 
    print ('\nPacking molecules...')
    load = raw_input('\nif packmol job is done, proceed to wrting topology files? (y) ')
    while load != 'y' and len(load)>0:
        print('enter y when packmol is done')
        load = raw_input('\nif packmol job is done, proceed to wrting topology files? (y) ')
    if load  == 'y':

#    for pdb in mixturePdb:
#        while not os.path.exists(pdb):  
#            time.sleep(1)  

        print('\nWriting input file to load force field in tleap')
        loadffFile,topFile,crdFile,pdbFile = loadFF.loadFF(watermodel,mixturePdb,PAHLib,ff)
        os.system('\ntleap -s -f {} > loadFF.out'.format(loadffFile))
    
    for i,top in enumerate(topFile):
        while not os.path.exists(top) and not os.path.exists(crdFile[i]):
           time.sleep(1)
    print "\nMaking sub-directories and move top and crd file to sub-directories"
    print "\nWriting sim.py files"
    for charge in f:
        for w in wtFraction:
            folder = 'f{}/w{}'.format(charge,w)
            if not os.path.exists(folder): 
                os.makedirs(folder)
            top = glob.glob('AH*f{}*w{}.parm7'.format(charge,w))
            crd = glob.glob('AH*f{}*w{}.crd'.format(charge,w))
	    pdb = glob.glob('AH*f{}*w{}.pdb'.format(charge,w))
            shutil.copy(top[0],folder+'/{}'.format(top[0]))
            shutil.copy(crd[0],folder+'/{}'.format(crd[0]))
	    shutil.copy(pdb[0],folder+'/{}'.format(pdb[0]))
            for file in pdb:
                shutil.copy(file,folder+'/{}'.format(file))
            for temp in T:
                sim = writeOpenMMinput(temp,top[0],crd[0],x,y,z,anisoP,np,N)
                job = writeJobfile(temp,charge,w,N)
                shutil.move(job,folder+'/{}'.format(job))
                shutil.move(sim,folder+'/{}'.format(sim))
    #remainingfiles = glob.glob('*')
    #if not os.path.exists('build'):
    #            os.makedirs('build')
    #for i in remainingfiles:
	#if os.path.isfile(i):
	#	shutil.move(i,'build')    
