#!/usr/bin/python3.4
import sys, mmap, re
from numpy import any
import numpy as np

def readOrbitals(infile):
    """ This function obtains a nwchem log-file (infile) and extracts the MO vectors from it (works only, 
    if in nwchem a odft-calculation was done). For this function, the patch 'movecs' is important for 
    working properly.
    """
    #load file and make it an mmap-object
    files=open(infile, "r")
    inp=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
    files.close
    #search for the alpha-part of MOvects
    atemp=re.findall(\
        b"(?<=DFT Final Alpha Molecular Orbital Analysis\n )[\w.=\+\- \n',^\"\d]+(?=DFT Final Beta)",
        inp, re.M)[-1]
    aMOvect=atemp.decode("utf-8").strip().split("Vector")
    anbf=len(aMOvect)-1 #because the first element is not an orbital vector
    anum,acoeff=getOrbitals(aMOvect[1:])
    #now, get the sorting and the first row to be printed
    aoccupation=getOcc(aMOvect[1:])
    aenergies=getEn(aMOvect[1:])
    
    # repeat for beta-porbitals
    btemp=re.findall(b"(?<=DFT Final Beta Molecular Orbital Analysis\n )[\d\w .=\+\- \n',^\"]+(?=\n\n)", inp, re.M)[-1]
    bMOvect=btemp.decode("utf-8").strip().split("Vector")
    bnbf=len(bMOvect)-1 
    bnum,bcoeff=getOrbitals(bMOvect[1:])
    boccupation=getOcc(bMOvect[1:])
    benergies=getEn(bMOvect[1:])
    
    # put other quantities in common vectors for returning
    occupation=[aoccupation, boccupation]
    energies=[aenergies, benergies]
    num=[anum,bnum]
    coeff=[acoeff,bcoeff]
    return num,coeff, occupation, energies

def getOrbitals( MOvect ):
    elements=[]
    NR=[]
    COEFF=[]
    for mo in range(len(MOvect)):
        currMOv=re.split('\n|            ', MOvect[mo].strip().split("---------------\n")[-1])
        orbital_nr=[]
        orbital_coeff=[]
        for i in range(len(currMOv)):
            elements.append(currMOv[i].split())
            orbital_nr.append(int(elements[i][0]))
            orbital_coeff.append(float(elements[i][1]))
        #resort elemnts by index
        NR.append(orbital_nr)
        COEFF.append(orbital_coeff)
    return NR,COEFF

def getEn(MOvect):
    NumOrbits=len(MOvect)
    Energies=np.zeros(NumOrbits)
    for ind in range(NumOrbits):
        currMOv=re.findall(r"(?<=E\=)[\d \.\+\-D]+", MOvect[ind])
        assert len(currMOv)==1, "There was an error reading occupation vector."
        Energies[ind]=float(currMOv[0].replace("D","e"))
    return Energies
    
def getOcc(MOvect):
    NumOrbits=len(MOvect)
    occupation=np.zeros(NumOrbits, dtype=int)
    for ind in range(NumOrbits):
        currMOv=re.findall(r"(?<=Occ\=)[\d]", MOvect[ind])
        assert len(currMOv)==1, "There was an error reading occupation vector."
        occupation[ind]=int(currMOv[0])
    return occupation
    
def print2molden(MO,coeff,occup, energy):
    #num,coeff, occupation, energies
    print("[MO]")
    for spin in [0,1]:
        for i in range(len(occup[spin])):
            print("Sym=c1") 
            print("Ene=%f"%energy[spin][i])
            if spin==0:
                print("Spin=Alpha")
            else:
                print("Spin=Beta")
            print("Occup=%d"%occup[spin][i])
            sort=np.argsort(MO[spin][i])
            for j in range(len(MO[spin][i])):
                print("%3d  %.7g"%(MO[spin][i][sort[j]],coeff[spin][i][sort[j]]))

def main(argv=None):
    # evaluate input-arguments
    assert len(argv)==1, "three input files expected."
    infile=argv[0] # input file
    
    # read first quantities from input-files
    
    num,coeff, occupation, energies=readOrbitals(infile)
    print2molden(num,coeff, occupation, energies)

if __name__ == "__main__":
   main(sys.argv[1:])

version=0.0

