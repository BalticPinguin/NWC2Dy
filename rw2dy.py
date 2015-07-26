#!/usr/bin/python
import re, mmap
import numpy as np

def getCoefficients( MOvect ):
   """ This function extracts the coefficients for molecular orbitals out-of log-files of
       a modified nwchem-version. The modification needed is: full print of MO-coefficients

       **PARAMETERS**
       MOvect    a vector containing the MO-coefficients as its elements;
                 each element contains a string for one MO.
                 To get the correct format, see function printOrbitals.
       **RETURNS**
       coeff     a matrix whose colums belong to MO vectors and rows contain coefficiens of
                 atomic orbitals.
   """
   for ind in range(len(MOvect)):
      currMOv=re.split('\n|            ', MOvect[ind].strip().split("---------------\n")[-1])
      elements=[]
      orbital_nr=[]
      for i in range(len(currMOv)):
         elements.append(currMOv[i].split())
         orbital_nr.append(int(elements[i][0]))

      #resort elemnts by index
      index=np.argsort(orbital_nr)
      if ind==0:
         coeff=np.zeros(( len(MOvect),len(elements) ))
      #fill elements into matrix
      for i in range(len(elements)):
         coeff[ind][i]=float(elements[index[i]][1])
   return coeff

def getFile(purpose):
   """ This function can be used to ask for additonal file names that are initially not given (e.g. due to
       extra options). 
       **PARAMETERS**
       purpose   a string containing a desciption of the file. It should tell the user,
                 which file has to be specified. (so make it descriptive ore even unique!)
       **RETURNS**
       the given name of the file

       **TODO**
       A good amendment to this function would be:
        a) to show availible files (ls)
        b) to make a test, whether the file exists and is readable
        c) probably even to make auto-verfollstaendigung; if possible in python??
   """
   infile = raw_input("Please type the file for %s :" %(purpose))
   return infile

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
      n=np.floor(pow(3*(shells[NA+1]-shells[NA]-1),1./3))+1
      assert not np.isnan(n), "The specification of atomic shells went wrong."
      n=int(n) #estimate number of shells.
      for a in range(int(shells[NA]), int(shells[NA])+n):
         orbitals=OrbitalNames(a-shells[NA]+1)
         j=0
         i=0
         while  i-j < min(len(orbitals), shells[NA+1]-electrons-shells[NA] ):
            #test, if orbitals coincide with those in  atom
            #print "start",i, electrons
            if orbitals[i][0][1:] !=atom[i-j+electrons][3]:
               #if they are not equal, take next orbital and compare with that one printed
               j+=1
            else:
               MO.append(orbitals[i][0])
               n=orbitals[i][1]
               l=orbitals[i][2]
               m=orbitals[i][3]
               neededSort.append(int(str(NA+1)+str(l)+str(m+l)+str(n))) #durty but works...
            i+=1
         electrons+=len(orbitals)-j
   aoe=np.argsort(neededSort)
   sortMO=[]
   for i in range(len(aoe)):
      sortMO.append("%s%s    %4s"%(atom[aoe[i]][1],atom[aoe[i]][2],MO[aoe[i]]))
   return sortMO, aoe

def getOcc(MOvect):
   NumOrbits=len(MOvect)
   occupation=np.zeros(NumOrbits, dtype=int)
   for ind in range(NumOrbits):
      currMOv=re.findall(r"(?<=Occ\=)[\d \.\+\-D]", MOvect[ind])
      assert len(currMOv)==1, "There was an error reading occupation vector."
      occupation[ind]=int(currMOv[0])
   return occupation

def OrbitalNames(n):
   """ wrong, if more than 21 shells have to be filled (after 8s, there will be errors)
       THIS FUNCTION HAS TO BE CHANGED MORE FUNDAMENTALLY TO BE MORE FLEXIBLE
       FOR DIFFERENT BASIS SETS.
       MOREOVER, IT SEEM TO CONTAIN CONTRADITORY INFORMATION ON THE MAIN QUANTUM NUMBER!
   """
   names=[]
   n=int(n)
   if n>=1:
      names.append([str(n)+"s",n ,0, 0])
   if n>=7:
      names.append([str(n-2)+"g4-",(n-1),4,-4])
      names.append([str(n-2)+"g3-",(n-1),4,-3])
      names.append([str(n-2)+"g2-",(n-1),4,-2])
      names.append([str(n-2)+"g1-",(n-1),4,-1])
      names.append([str(n-2)+"g0",(n-1),4, 0])
      names.append([str(n-2)+"g1+",(n-1),4, 1])
      names.append([str(n-2)+"g2+",(n-1),4, 2])
      names.append([str(n-2)+"g3+",(n-1),4, 3])
      names.append([str(n-2)+"g4+",(n-1),4, 4])
   if n>=6:
      names.append([str(n-2)+"f3-",(n-1),3,-3])
      names.append([str(n-2)+"f2-",(n-1),3,-2])
      names.append([str(n-2)+"f1-",(n-1),3,-1])
      names.append([str(n-2)+"f0",(n-1),3, 0])
      names.append([str(n-2)+"f1+",(n-1),3, 1])
      names.append([str(n-2)+"f2+",(n-1),3, 2])
      names.append([str(n-2)+"f3+",(n-1),3, 3])
   if n>=4: 
      names.append([str(n-1)+"d2-",(n-2),2,-2])
      names.append([str(n-1)+"d1-",(n-2),2,-1])
      names.append([str(n-1)+"d0",(n-2),2, 0])
      names.append([str(n-1)+"d1+",(n-2),2, 1])
      names.append([str(n-1)+"d2+",(n-2),2, 2])
   if n>=2:
      names.append([str(n)+"px",n,1,-1])
      names.append([str(n)+"py",n,1,0])
      names.append([str(n)+"pz",n,1,1])
   return names 

def printCI(outfile, CIcoeff, transition,sort, noocc, nofree, states,occ):
   output=open(outfile, 'a')
   statestr=""
   for i in range(noocc):
      if occ[i]==1:
         statestr+="u"
      else:
         statestr+=str(occ[i])
   for i in range(nofree):
      if occ[i]==1:
         statestr+="d"
      else:
         statestr+=str(occ[i+noocc])
   for i in range(states):
      output.write("\n\nSTATE=%d \n" %(i+1))
      output.write("Det             Occupation  ")
      for nomatter in range(noocc+nofree-10):
         output.write(" ")
      output.write("    Coef                      Weight\n")
      for j in range(noocc*nofree):
         #print first (counting) number
         output.write("  %3d         "%(j+1))
         #print the occupation of the state
         sortstate=""
         for k in range(noocc+nofree):
            #print transition[j]
            if k==transition[j][0]-1:
               if occ[k]==1:
                  sortstate+="0"
               elif occ[k]==2:
                  sortstate+="d"
            elif k==transition[j][1]-1:
               if occ[k]==0:
                  sortstate+="u"
               elif occ[k]==1:
                  sortstate+="2"
            else:
               sortstate+=statestr[k]
         output.write("%s"%sortstate)
         #print coefficient and weight
         output.write("     %16.10g         %16.10g\n"%( CIcoeff[i][j], CIcoeff[i][j]*CIcoeff[i][j]))
   output.close()

def printnorm(CIcoeff):
   nos=len(CIcoeff)
   trans=len(CIcoeff[0])
   for i in range(nos):
      norm=0
      for j in range(trans):
         norm+=CIcoeff[i][j]*CIcoeff[i][j]
      print norm

def printOrbitals(infile,outfile):
   output=open(outfile, 'a')
   files=open(infile, "r")
   inp=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
   files.close
   #temp=re.findall(r"(?<=ROHF Final Molecular Orbital Analysis\n )[\d\w .=\+\- \n',^\"]+(?=center of mass)", inp, re.M)[-1]
   temp=re.findall(r"(?<=DFT Final Molecular Orbital Analysis\n )[\d\w .=\+\- \n',^\"]+(?=center of mass)",  inp, re.M)[-1]
   MOvect=temp.strip().split("Vector")
   nbf=len(MOvect)-1 #because the first element is not an orbital vector
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
               output.write(u"  %0.10E " %(indMatrix[i][sort[j]])) #convert the repective element
            output.write("\n")
   output.close
   return sort

def readCI(CIfile):
   #open the file 
   cifi=open(CIfile, "r")
   #throw away first two lines
   cifi.readline()
   cifi.readline()
   #number of states
   nos=int(float(cifi.readline()))
   #number of occupied orbitals
   noocc=int(float(cifi.readline().strip().split()[0]))
   #total number of orbitals
   notot=int(float(cifi.readline().strip().split()[0]))
   nouno=notot-noocc #number of unoccupied orbitals
   #throw away next 4 lines
   for i in range(4):
      cifi.readline()

   #now, start serious work:
   CIcoeff=np.zeros((nos,nouno*noocc))
   CItrans=np.zeros((nouno*noocc,2))
   for state in range(nos):
      cifi.readline()
      cifi.readline()
      for trans in range(nouno*noocc):
         CIcoeff[state][trans]=float(cifi.readline())
         #print " %2d   %.7E   %3d  %3d" %(state, CIcoeff[state][trans], trans/(nouno)+1, noocc+trans%(nouno)+1)
         if state==0:
            CItrans[trans]=[trans/(nouno), trans%(nouno)+noocc]
   cifi.close()
   return CIcoeff, CItrans, noocc, nouno, nos

def readCI2(infile):
   #open the file 
   files=open(infile, "r")
   inp=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
   files.close()
   
   roots=re.findall(r"(?<=Root )[ \d]+", inp)
   for i in range(len(roots)):
      roots[i]=int(roots[i])
   nos=np.max(roots)
   #cicoeff=re.findall(r"(?<=Dipole Oscillator Strength )[Occ. Virt. abE\-\+\d\n]+", inp)
   cicoeff=re.findall(r"(?<=Dipole Oscillator Strength )[Occ. Virt. abE\-\+\d\n alpha beta]+", inp)
   if len(cicoeff)>nos:
      cicoeff=cicoeff[-nos:] #this gives only last ci-vector
   elif len(cicoeff)<nos:
      assert 1==2, "an error occured. Not all roots found."
   for i in range(nos):
      CI=cicoeff[i].strip().split("\n")[2:] #split into lines and through first and last line away
      if i==0:
         if "Occ." in CI[-1]:
            noorb=len(CI)
         else:
            noorb=len(CI)-1
         CIcoeff=np.zeros((nos,noorb))
         CItrans=np.zeros((noorb,2))
      for j in range(noorb):
         transition=CI[j].split()
         CIcoeff[i][j]=transition[-1]
         if i==0:
            if len(transition)==8:
               CItrans[j]=[transition[1],transition[5]]
            elif len(transition)==10:
               CItrans[j]=[transition[1],transition[6]]
            else:
               assert 1==2, "an error occured. Structure of CI-vectors unknown: \n %s\n" %transition
   #print len(CItrans)
   noocc=int(np.max(CItrans[:].T[0]))
   nouno=int(np.max(CItrans[:].T[1]))-noocc
   return CIcoeff, CItrans, noocc, nouno, nos

def readOrbitals(infile):
   files=open(infile, "r")
   inp=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
   files.close
   #temp=re.findall(r"(?<=ROHF Final Molecular Orbital Analysis\n )[\d\w .=\+\- \n',^\"]+(?=center of mass)", inp, re.M)[-1]
   temp=re.findall(r"(?<=DFT Final Molecular Orbital Analysis\n )[\d\w .=\+\- \n',^\"]+(?=center of mass)", inp, re.M)
   assert temp!=[], "Either no DFT MO analysis is taken into account or the calculation terminated with an error."
   temp=re.findall(r"(?<=DFT Final Molecular Orbital Analysis\n )[\d\w .=\+\- \n',^\"]+(?=center of mass)", inp, re.M)[-1]
   MOvect=temp.strip().split("Vector")
   nbf=len(MOvect)-1 #because the first element is not an orbital vector
   #now, fet the sorting and the first row to be printed
   MOs, sort=getOrbitals( MOvect[1] )
   occupation=getOcc(MOvect[1:])
   return MOs, sort, occupation

def rOverlap(infile, sort):
   files=open(infile, "r")
   inp=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
   files.close

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
   for i in range(len(indices)):
      overlap[indices[i][1]-1][indices[i][0]-1]=values[i] #make symmetric matrix
   #  -> overlap /overlap 2 differ (by small numbers)
   return overlap,dim

def rwenergy(output, infile):
   files=open(infile, "r")
   inp=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
   files.close()
   energy=float(re.findall(r"(?<=Total SCF energy =)[\d \-\.]+", inp)[-1])
   output.write("%15.10g\n"%energy)
   exen=re.findall(r"(?<=Root )[\w \d\.\-]+ ",inp)
   roots=[]
   for i in range(len(exen)):
      if "." in exen[i]:
         roots.append(float(exen[i].strip().split()[0]))
   nos=int(np.max(roots))
   if "a.u." in exen[-1]:
      last=-1
   else:
      last=-2
   exen=exen[-nos+last:last] 
   for i in range(len(exen)):
      output.write("%15.10g\n"%(float(exen[i].strip().split()[-4])+energy))

def writePreamble(outfile, noocc, nouno, nbf):
   output=open(outfile, 'w')
   output.write("MOLCAS\n0\n\n")
   output.write("METHOD\n1\n\n")
   output.write("FNSPIN\n1\n\n")
   output.write("INSPIN\n1\n\n")
   output.write("FMULT\n2\n\n")
   output.write("IMULT\n1\n\n")
   output.write("NACTIVE\n%d\n\n"%(noocc+nouno))
   output.write("NFROZEN\n0\n\n")
   output.write("FNACTEL\n%d\n\n"%(2*noocc-1))
   output.write("INACTEL\n%d\n\n"%(2*noocc))
   output.write("NBASF\n%d\n\n"%(nbf))
   output.write("SFPRINT\n3\n\n")
   output.write("SOCPRINT\n0\n\n")
   output.write("FINALWF\n%d  %d\n\n"%(noocc*nouno,noocc*nouno))
   output.write("INITIALWF\n%d  %d\n\n"%(noocc*nouno,noocc*nouno))
   #output.write("INITIALWF\n1  1\n\n")
   output.close()

def writePreamble2(outfile, noocc,nouno, nbf, occi,occf):
   output=open(outfile, 'w')
   output.write("MOLCAS\n0\n\n")
   output.write("METHOD\n1\n\n")
   output.write("FNSPIN\n1\n\n")
   output.write("INSPIN\n1\n\n")
   output.write("FMULT\n2\n\n")
   output.write("IMULT\n1\n\n")
   output.write("NACTIVE\n%d\n\n"%(noocc+nouno))
   output.write("NFROZEN\n0\n\n")
   output.write("FNACTEL\n%d\n\n"%(occf))
   output.write("INACTEL\n%d\n\n"%(occi))
   output.write("NBASF\n%d\n\n"%(nbf))
   output.write("SFPRINT\n3\n\n")
   output.write("SOCPRINT\n0\n\n")
   output.write("FINALWF\n%d  %d\n\n"%(noocc*nouno,noocc*nouno))
   #output.write("INITIALWF\n%d  %d\n\n"%(noocc*nouno,noocc*nouno))
   output.write("INITIALWF\n1  1\n\n")
   output.close()

def wOverlap(overlap, dim,sort, outfile):
   output=open(outfile, 'a')
   output.write("\nATOMOVERLAP\n")
   #print matrix in required format:
   for n in range(0,dim,5):
         output.write(u"\n")
         for i in range(n,min(n+5,dim)):
            output.write(u"                %2i" %(i+1))
         output.write("\n")
         for j in range(dim):
            output.write(u"  %2i  "%(j+1))
            for i in range(n,min(n+5,dim)):
               #output.write(u"  %0.10E " %(overlap[i][j])) #convert the repective element
               output.write(u"  %0.10E " %(overlap[sort[i]][sort[j]])) #convert the repective element
            output.write("\n")
   output.close

