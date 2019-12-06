#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:06:15 2019

@author: nvthaomy
"""

def loadFF(watermodel,mixturePdb,lib,ff):
    """write tleap input file to load forcefield and generate .parm7 and .crd
        mixture Pdb is list of all mixture pdb files
	lib: name of tleap library for PAH monomers"""
    with open('loadFF.in','w') as load:
        load.write('source leaprc.{}'.format(ff))
        load.write('\nsource leaprc.water.{}'.format(watermodel))
        for i in lib:
            load.write('\nloadOFF {}.lib\n'.format(i))
        topFile = []
        crdFile = []
        pdbFile = []
        for pdb in mixturePdb:
            topname = pdb[:pdb.index('w')]+str(ff)+'_'+pdb[pdb.index('w0'):pdb.index('pdb')]+'parm7'
            crdname = pdb[:pdb.index('pdb')]+'crd'
            pdbname = topname.split('.parm7')[0]+'.pdb'
            topFile.append(topname)
            crdFile.append(crdname)
            pdbFile.append(pdbname)
            load.write('\n\nx=loadpdb {}'.format(pdb))
            load.write('\naddions x Cl- 0')
            load.write('\nsetbox x vdw 1')
            load.write('\nsaveamberparm x {} {}'.format(topname,crdname))
	    load.write('\nsavepdb x {}'.format(pdbname))
        load.write('\nquit') 
    return 'loadFF.in',topFile,crdFile,pdbFile
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p",nargs='+', required=True,
                        help="list of mixture pdb")
    parser.add_argument("-w","--watermodel", default = 'opc',
                        help="Water model for simulation (opc,tip3p,spce), default = opc")
    parser.add_argument("-l", nargs='+', defaul = ['PAH_deprot.lib','PAH_prot.lib'],
			help="tleap libraries for PAH monomers: PAH, PAH_avg, PAH1, etc.")
    parser.add_argument("-ff", default = 'ff99',
                        help="amber forcefield that compatiple with atom types in libraries")
    args = parser.parse_args() 
    watermodel = args.watermodel
    ff = args.ff
    mixturePdb = args.p
    watermodel = args.w
    loadFF(watermodel,mixturePdb,args.l,ff)
