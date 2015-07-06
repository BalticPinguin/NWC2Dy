#!/usr/bin/python
import re, mmap, sys, decimal
import numpy as np

def readPreamble(preamble, outfile):
   output=open(outfile, 'w')
   with open(preamble, "r") as pream:
       for line in pream.readlines():
         output.write("%s\n" %line)
   output.close

def rwOverlap(infile, outfile):
   output=open(outfile, 'a')
   files=open(infile, "r")
   inp=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
   files.close

   output.write("ATOMOVERLAP\n")
   temp=re.findall(r"(?<=\n Begin overlap 1-e integrals\n ====================================\n)[\d.\w \+\- \n]+", inp)[-1]
   lines=temp.strip().split("\n")
   indices=[]
   values=[]
   for i in range(len(lines)):
      foo=re.findall(r"(?<=  )[\d]+",lines[i])        #gives 2 numbers per element: i,j
      bar=[int(foo[0]) ,int(foo[1])] #select first two elements each
      indices.append(bar)  # and put into index-vector
      foo=re.findall(r" [- ]0.[\d\.E\-\+]+",lines[i])  #gives the overlap each
      bar=float(foo[0])
      values.append(bar)
   #now, write into matrix:
   dim=np.max(indices)
   overlap=np.zeros((dim,dim))
   overlap2=np.zeros((dim,dim))
   for i in range(len(indices)):
      overlap2[indices[i][0]-1][indices[i][1]-1]=values[i]
      overlap[indices[i][1]-1][indices[i][0]-1]=values[i] #make symmetric matrix
   #  -> overlap /overlap 2 differ (by small numbers)

   #print matrix in required format:
   for n in range(0,dim,5):
         output.write(u"\n")
         for i in range(n,min(n+5,dim)):
            output.write(u"                %2i" %(i+1))
         output.write("\n")
         for j in range(dim):
            output.write(u"  %2i  "%(j+1))
            for i in range(n,min(n+5,dim)):
               output.write(u"  %0.10E " %(overlap[i][j])) #convert the repective element
            output.write("\n")
   output.close

def printOrbitals(infile,outfile):
   output=open(outfile, 'a')
   files=open(infile, "r")
   inp=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
   files.close
   temp=re.findall(r"(?<=ROHF Final Molecular Orbital Analysis\n )[\d\w .=\+\- \n',^\"]+(?=center of mass)", inp, re.M)[-1]
   MOvect=temp.strip().split("Vector")
   nbf=len(MOvect)-1 #because the first element is not an orbital vector
   #now, fet the sorting and the first row to be printed
   MOs, sort=getOrbitals( MOvect[1] )
   #obtain the coefficients
   indMatrix=getCoefficients(MOvect[1:])

   #print matrix in required format:
   for n in range(0,nbf,5):
         output.write(u"\n\n       Orbital")
         "                %i"
         for i in range(n,min(n+5,nbf)):
            output.write(u"                %i" %(i+1))
         output.write("\n\n")
         for j in range(len(MOs)):
            output.write(u"  %2i    %s"%(j+1, MOs[j]))
            for i in range(n,min(n+5,nbf)):
               output.write(u"  %0.10E " %(indMatrix[sort[i]][j])) #convert the repective element
            output.write("\n")
   output.close

def OrbitalNames(n):
   """ wrong, if more than 21 shells have to be filled (after 8s, there will be errors)
   """
   names=[]
   n=int(n)
   if n>=1:
      names.append([str(n)+"s  ",n ,0, 0])
   if n>=6:
      names.append([str(n-2)+"f3-",(n-1),3,-3])
      names.append([str(n-2)+"f2-",(n-1),3,-2])
      names.append([str(n-2)+"f1-",(n-1),3,-1])
      names.append([str(n-2)+"f0 ",(n-1),3, 0])
      names.append([str(n-2)+"f1+",(n-1),3, 1])
      names.append([str(n-2)+"f2+",(n-1),3, 2])
      names.append([str(n-2)+"f3+",(n-1),3, 3])
   if n>=4:
      names.append([str(n-1)+"d2-",(n-2),2,-2])
      names.append([str(n-1)+"d1-",(n-2),2,-1])
      names.append([str(n-1)+"d0 ",(n-2),2, 0])
      names.append([str(n-1)+"d1+",(n-2),2, 1])
      names.append([str(n-1)+"d2+",(n-2),2, 2])
   if n>=2:
      names.append([str(n)+"px ",n,1,-1])
      names.append([str(n)+"py ",n,1,0])
      names.append([str(n)+"pz ",n,1,1])
   return names 

def getCoefficients( MOvect ):
   for ind in range(len(MOvect)):
      currMOv=re.split('\n|            ', MOvect[ind].strip().split("---------------\n")[-1])
      elements=[]
      orbital_nr=[]
      for i in range(len(currMOv)):
         elements.append(currMOv[i].split())
      if ind==0:
         coeff=np.zeros(( len(MOvect),len(elements) ))

      #resort elemnts by index
      for i in range(len(elements)):
         coeff[ind][i]=float(elements[i][1])
   return coeff

def getOrbitals( MOvect ):
   currMOv=re.split('\n|            ', MOvect.strip().split("---------------\n")[-1])
   elements=[]
   orbital_nr=[]
   for i in range(len(currMOv)):
      elements.append(currMOv[i].split())
      orbital_nr.append(int(elements[-1][0]))

   #resort elemnts by index
   index=np.argsort(orbital_nr)
   NumAtom=int(elements[index[-1]][2])
   shells=np.zeros(NumAtom+1)
   atom=[] #contains the name of an atom and the atomic orbital
   temp=0
   shells[0]=0
   for i in range(len(index)):
      atom.append([elements[index[i]][0], elements[index[i]][2], elements[index[i]][3],elements[index[i]][4]])
      if temp+1!=int(atom[-1][1]):
         temp+=1         #purpose of this: keep track, where next atom starts
         shells[temp]=i
   shells[-1]=len(index)
   MO=[]
   neededSort=[]
   for NA in range(NumAtom):
      electrons=0
      for a in range(int(shells[NA]), int(shells[NA+1])):
         orbitals=OrbitalNames(a-shells[NA]+1)
         for i in range(int(min(len(orbitals), shells[NA+1]-electrons-shells[NA] ) )):
            MO.append(orbitals[i][0])
            n=orbitals[i][1]
            l=orbitals[i][2]
            m=orbitals[i][3]
            neededSort.append(int(str(NA+1)+str(l)+str(m+l)+str(n))) #durty but works...
         electrons+=len(orbitals)
   aoe=np.argsort(neededSort)
   sortMO=[]
   totalsort=[]
   for i in range(len(aoe)):
      sortMO.append("%s%s    %s"%(atom[aoe[i]][1],atom[aoe[i]][2],MO[aoe[i]]))
      totalsort.append(index[aoe[i]])
   return sortMO, totalsort

def main(argv=None):
   assert len(argv)==3, "three input files expected."
   preamble=argv[0]
   infile=argv[1]
   outfile=argv[2]
   readPreamble(preamble, outfile)
   rwOverlap(infile,outfile)
   printOrbitals(infile,outfile)

if __name__ == "__main__":
   main(sys.argv[1:])

version=0.1
