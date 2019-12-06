#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 17:56:07 2019

@author: nvthaomy
"""
import random
import numpy as np
import math
import os

""" write out tleap input file to build a single PAH chain with specified
    tacticity and various fraction of ionization
    input: N: degree of polymerization
           f: fraction of ionization
           Pm: fraction of meso diads"""
class Molecule:
    def __init__(self,N,f,Pm):
        self.DOP = N
        self.charge = f #a list of ionization fraction e.g. [0,0.2,0.6,1]
        self.Pm = Pm
    def get_indices(self,num1,num2,N):
        """get indices of ionized monomers (for f<=0.5) or ionizedonated monomer
           (for f>0.5) in a chain
           N: degree of polymerization
           num2: number of ionized monomers (for f<=0.5) or ionizedonated monomer (for f>0.5)
           num1: list of upper and lower bounds of number of ionizedonated monomers if f  <= 0.5 or 
                number of ionized monomers if f > 0.5
           output: list of indices"""
        i=[] #list of indices of ionized monomers
        if len(num1) > 1: #for case when num is not an integer or divisible by 5
            i.append(random.randint(0,min(num1)))
            print "place the first neutral (or ionized if f>0.5) monomer at position %i" %(i[0]+1)
            counter = 0 
            old_j = 0
            for j in range(1,int(num2)):
                spacing = num1[random.randint(0,len(num1)-1)] #draw random spacing from the list of upper and lower bounds of number of ionizedonated(ionized) monomers 
                i.append(i[j-1] + spacing + 1)
                while i[j] > N-1:
                    i[j] -= 1        
                    if j != old_j: #only count if the next monomer index exceed the monomer chain length
                        counter += 1 
                        if counter > 1:
                            print "\nWarning: Index of monomer is out of range!"
                    old_j = j 
        else:
            i.append(0)
            counter = 0
            old_j = 0
            print "place the first neutral/ionized monomer at position %i" %(i[0]+1)
            for j in range(1,int(num2)):
                i.append(i[j-1] + num1[0] + 1)
                while i[j] > N-1:
                        i[j] -= 1        
                        if j != old_j: #only count if the next monomer index exceed the monomer chain length
                            counter += 1 
                            if counter > 1:
                                print "\nWarning: Index of monomer is out of range!"
                        old_j = j 
        i = [int(x) for x in i] #convert indices to integer
        return i
            
    def uniform_charge_pattern(self,charge_frac):
        print "\nGenerating evenly distributed charge pattern for f = %3.2f" %charge_frac 
        N=self.DOP
        f0= N*['d']
        f1= N*['p']
        
        if charge_frac == 0:
            charge_pattern = f0
        if charge_frac == 1:
            charge_pattern = f1
        if charge_frac != 0 and charge_frac != 1:
            print "\nDeionizedonation fraction is not 0 or 1"
            
            if charge_frac <= 0.5:
                num_neutral = (charge_frac * N)
                num_ionized = [(N/num_neutral - 1)] 
                print "\nNumber of ionizedonated monomers between each ionized monomer is: %3.3f" %num_ionized[0]
                if not abs(num_ionized[0]-round(num_ionized[0])) < 0.0001: #not an integer
                    if num_ionized[0]*10 % 5 == 0: #screen out charge fraction that result in num_ionized= x.5
                        num_ionized = [math.ceil(num_ionized[0])]
                    else:
                        num_min = math.floor(num_ionized[0])
                        num_max = math.ceil(num_ionized[0])
                        mean = num_ionized[0] #average number of ionizedonated moners in between two ionized ones
                        num_ionized = [num_min]
                        calc_mean = sum(num_ionized)/len(num_ionized)
                        
                        while abs(calc_mean - mean) > 0.001:
                            if calc_mean > mean:
                                num_ionized.append(num_min)
                                calc_mean = sum(num_ionized)/len(num_ionized)
                            else:
                                num_ionized.append(num_max)
                                calc_mean = sum(num_ionized)/len(num_ionized)
                        print calc_mean
                        print ('\nNunmber of ionizedonated monomers will be drawn from this list:')
                        print num_ionized
                
                i_neutral = self.get_indices(num_ionized,num_neutral,N)
                print "\n Indices of ionized monomer for f = %3.2f are " %charge_frac
                print i_neutral
                if abs(len(i_neutral)-num_neutral) > 0.001:
                    print "Number of ionized monomers does not match the value of charge fraction!\n"
                charge_pattern = f0
                for i in i_neutral:
                    charge_pattern[i] = 'p'
            
            if charge_frac > 0.5:
                num_ionized = ((1 - charge_frac) * N)
                num_neutral = [(N/num_ionized - 1)]
                print "\nNumber of ionized monomers between each ionizedonated monomer is: %3.3f" %num_neutral[0]
                if not abs(num_neutral[0]-round(num_neutral[0])) < 0.0001: #not an integer
                    if num_neutral[0]*10 % 5 == 0: #screen out charge fraction that result in num_neutral= x.5
                        num_neutral = [math.ceil(num_neutral[0])]
                    else:
                        num_min = math.floor(num_neutral[0])
                        num_max = math.ceil(num_neutral[0])
                        mean = num_neutral[0] #average number of ionized moners in between two ionizedonated ones
                        num_neutral = [num_min]
                        calc_mean = sum(num_neutral)/len(num_neutral)
                        
                        while abs(calc_mean - mean) > 0.001:
                            if calc_mean > mean:
                                num_neutral.append(num_min)
                                calc_mean = sum(num_neutral)/len(num_neutral)
                            else:
                                num_neutral.append(num_max)
                                calc_mean = sum(num_neutral)/len(num_neutral)
                        print calc_mean
                        print ('\nNunmber of ionized monomers will be drawn from this list:')
                        print num_neutral
                
                i_ionized = self.get_indices(num_neutral,num_ionized,N)
                print "\nIndices of ionizedonated monomer for f = %3.2f are " %charge_frac
                print i_ionized
                if abs(len(i_ionized)-num_ionized) > 0.001:
                    print "Number of ionizedonated monomers does not match the value of charge fraction!\n"
                charge_pattern = f1
                for i in i_ionized:
                    charge_pattern[i] = 'd' 
        return charge_pattern
    
    def charge_pattern(self,charge_frac,pattern):
        """Evaluate if charge pattern is random or evenly distributed
        and enerate charge pattern"""
        if pattern == 'random':
            print "\nGenerating randomly distributed charge pattern for f = %3.2f" %charge_frac 
            N=self.DOP
            f0= N*['d']
            f1= N*['p']
            num_neutral = N*charge_frac
            if charge_frac == 0:
                charge_pattern = f0
            if charge_frac == 1:
                charge_pattern = f1
            else:
                charge_pattern = f0
                i_neutral = random.sample(range(N),int(num_neutral))
                print i_neutral
                for i in i_neutral:
                    charge_pattern[i] = 'p' 
            print charge_pattern
        else:
            charge_pattern = self.uniform_charge_pattern(charge_frac)
        return charge_pattern           
                
        
    def join_tact_charge(self,charge_frac,tacticity,pattern):
        """Append tacticity with charge pattern"""
        tact_charge=[]
        charge_pattern = self.charge_pattern(charge_frac,pattern)
        for i in range(0,self.DOP):
            a = tacticity[i] + charge_pattern[i]
            tact_charge.append(a)
        tact_charge[0]=tact_charge[0].upper()
        tact_charge[-1]=tact_charge[-1].upper()
        return tact_charge
        
    def tacticity(self,Pm):
            N=self.DOP
            dyad=[]
            tact=[]
            #follow Bernoullian statistics (common with free radical polymerization)
            #stereo of the next monomer is independent of the stereochemistry of growing chain
            for i in range(0,N):
                rand = np.random.random()
                if rand <= Pm:
                        dyad.append('m')
                else:
                        dyad.append('r')
            print('Diad sequence from meso diad fraction of {}:\n{}'.format(Pm,dyad))
            tact.append('NH') #head group is achiral
            tact.append('u') #arbitrarily pick the second monomer to be "up"
            for i in range(1,len(dyad)-2):
                if dyad[i] == 'm':
                    tact.append(tact[i])
                else:
                    if tact[i] == 'u':
                        tact.append('d')
                    else:
                        tact.append('u')
            tact.append('NT')
            return tact
def buildPAH(f,N,Pm,pattern):
    with open("build_AH"+str(N)+".log", 'w') as log:     
        PAH=Molecule(N,f,Pm) #f is a list of ionization fraction
        tact_charge_matrix=[] #each row is the charge pattern of each ionizedonation fraction
        log.write("\nCalculating number of ionized monomers and"\
                  "modifying input fraction of ionization:")        
        for index,charge_frac in enumerate(f):
            if charge_frac !=0:
                num_neutral = float(charge_frac * N)
                if not num_neutral.is_integer(): 
                    log.write("\nNeed to modify f = %3.2f to " %charge_frac)
                    num_neutral = round(num_neutral)
                    charge_frac = num_neutral/N #new charge fraction
                    log.write("f = %3.2f" %charge_frac)
                    f[index] = charge_frac
        log.write("\nNew charge vector is :"%f)
        rounded_charge = [round(charge,2) for charge in f]
        log.write('{}'.format(rounded_charge))
        log.write('\nBuilding PAH of with N = {}, f = {}, charge pattern = {},'\
                  ' meso fraction  = {}'.format(N,rounded_charge,pattern,Pm))
        tacticity=PAH.tacticity(Pm)
        for index,charge_frac in enumerate(f):
            tact_charge = PAH.join_tact_charge(charge_frac,tacticity,pattern)
            tact_charge_matrix.append(tact_charge)
        
        log.write("\nWriting tleap input file to build a single polymer with different ionized fraction")
        file = open("build_AH"+str(N)+".in","w")
        file.write("source leaprc.gaff2\n")
        file.write("source leaprc.ff99\n") 
        file.write("loadOFF PAH_deprot.lib\n")
        file.write("loadOFF PAH_prot.lib\n")
        file.write("up =loadpdb PAH_prot1.pdb\n")
        file.write("ud =loadpdb PAH_deprot1.pdb\n")
        file.write("dp =loadpdb PAH_prot2.pdb\n")
        file.write("dd =loadpdb PAH_deprot2.pdb\n")
        file.write("set up head up.1.1\n")
        file.write("set up tail up.1.3\n")
        file.write("set ud head ud.1.1\n")
        file.write("set ud tail ud.1.3\n")
        file.write("set dp head dp.1.1\n")
        file.write("set dp tail dp.1.3\n")
        file.write("set dd head dd.1.1\n")
        file.write("set dd tail dd.1.3\n")
        file.write("\n")
        pdbList = []
        newCharge = []
        for index,charge_frac in enumerate(f):
            if abs(charge_frac - 0) < 10**(-3) or abs(charge_frac - 1) < 10**(-3):
                charge_frac = int(charge_frac)
            else:
                charge_frac = round(charge_frac,2)
            newCharge.append(charge_frac)
            sequence =' '.join(tact_charge_matrix[index])
            file.write('#f = {}\n'.format(round(charge_frac,2)))
            file.write("x = sequence{")
            file.write("{}".format(sequence))
            file.write("}\n")
            file.write("savepdb x AH{}_f{}.pdb\n".format(N,charge_frac))
            file.write("\n")  
            pdbList.append('AH{}_f{}.pdb'.format(N,charge_frac))
        file.write("quit")
        file.close()
        log.write("\nDone writing tleap input file")
    log.close()
    build_file = "build_AH"+str(N)+".in"
    return build_file,pdbList,newCharge
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",nargs='+', type=float, required=True,
                        help="list of ionization fractions e.g. 0 0.2 0.5 1")
    parser.add_argument("-N",type=int, required=True,
                        help="degree of polymerization")
    parser.add_argument("-Pm","--Pm",type=float,default = 1, 
                        help="Meso diad fraction, isotactic if = 1, syndiotactic if = 0")
    parser.add_argument("-r","--random", action = "store_true",
                        help="Random ionization, default pattern is random")
    parser.add_argument("-e","--evendist", action = "store_true",
                        help="Evenly distributed ionization")
    args = parser.parse_args()    
    f = args.f
    N = args.N
    Pm = args.Pm
    pattern = 'random' #random ionization
    if args.evendist:
        pattern = 'even'
    buildPAH(f,N,Pm,pattern)


     



