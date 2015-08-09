#!/usr/bin/python3
import re, mmap
import numpy as np

# changelog 0.3
# 1)
# changelog 0.2
# 1) converted to python3
# 2) got running code for unrestricted system

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

def getOrbitals2( MOvect ):
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

def getOrbitals( MOvect ):
   currMOv=re.split('\n|            ', MOvect.strip().split("---------------\n")[-1])
   elements=[]
   orbital_nr=[]
   for i in range(len(currMOv)):
      elements.append(currMOv[i].split())
      #print(elements[i])
      orbital_nr.append(int(elements[i][0]))
   index=np.argsort(orbital_nr)
   #resort elemnts by index
   MO=[]
   for i in range(len(index)):
      # this is a test. I don't know, if it will work...
      if len(elements[index[i]])>5:
         #print(elements[index[i]][4:])
         elements[index[i]][4]=elements[index[i]][4]+elements[index[i]][5]
         #print(elements[index[i]][4])
      MO.append([elements[index[i]][2], elements[index[i]][3],elements[index[i]][4]])
      #print(MO[i])
   #print("\n\n")
   NumAtom=int(MO[-1][0])
   sort=[]
   sortMO=[]
   # this is a test. I don't know, if it will work...
   control=False
   for foo in elements:
      if re.search("f", foo[4]) is not None:
         types=["s", "px", "py", "pz","d-2", "d-1", "d0", "d1", "d2","f-3", "f-2", "f-1", "f0", "f1", "f2","f3"]
         n=[0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3]
         control=True
         break
   if not control:
      for foo in elements:
         if re.search("d", foo[4]) is not None:
            types=["s", "px", "py", "pz","d-2", "d-1", "d0", "d1", "d2"]
            n=[0, 1, 1, 1, 2, 2, 2, 2, 2]
            control=True
            break
   if not control:
      types=["s", "px", "py", "pz"]
      n=[0, 1, 1, 1]
   while len(sort)<len(MO):
      #yes, here I do use an old name for a new object.
      for atom in range(NumAtom): 
         for x in range(len(types)):
            m=0
            for i in range(len(MO)):
               # this is not very efficient but ...
               if int(MO[i][0])-1!=atom:
                  continue
               if MO[i][2]==types[x]:
                  m+=1
                  sort.append(i)
                 # print(MO[i], n[x])
                  sortMO.append("%s%s    %i%s"%(MO[i][0], MO[i][1], n[x]+m, MO[i][2]))
                  #print(sort[-1], sortMO[-1], orbital_nr[index[i]])
            if len(sort)==len(elements): #all functions found
               break #stop the iteration (ready)
   #print(sortMO)
   #print(sort)
 #  for i in range(len(sort)):
 #     print("   ", sort[i], sortMO[i])
 #  print("\n\n")
   return sortMO, sort

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

def printCI(outfile, CIcoeff, transition,sort, noocc, nofree, states,occ, trans):
   def Trans(transition, state):
      newstate=""
      sortstate=Trans(transition[j], statestr)
      for i in range(state):
         if str(i+1) in transition:
           #ckeck, which one it is.
           occ[0][transition[2]-1]="1" #: --> work with occ[]?
         else:
            newstate+=state[i]
      return state

   output=open(outfile, 'a')
   statestr=""
   for i in range(noocc+nofree):
      if occ[0][i]==1 and occ[1][i]==1:
         statestr+="2"
      elif occ[0][i]==1:
         statestr+="u"
      elif occ[1][i]==1:
         statestr+="d"
      else:
         statestr+="0"

   for i in range(states):
      output.write("\n\nSTATE=%d \n" %(i+1))
      output.write("Det             Occupation  ")
      for nomatter in range(noocc+nofree-10):
         output.write(" ")
      output.write("    Coef                      Weight\n")
      for j in range(trans): #trans is number of transitions
         #print first (counting) number
         output.write("  %3d         "%(j+1))
         #print the occupation of the state
         #sortstate=Trans(transition[j], statestr, occ)
         if transition[j][0]!=0: #a transition from occ. alpha-state:
            if statestr[transition[j][0]-1]=="2":
               #######======= this is very slow.
               sortstate="".join("d" if i==transition[j][0]-1 else c for i,c in enumerate(statestr))
            elif occ[0][transition[j][0]-1]==1:
               sortstate="".join("0" if i==transition[j][0]-1 else c for i,c in enumerate(statestr))
            else:
               print(occ[0][transition[j][0]-1], transition[j][2])
         elif transition[j][2]!=0: #a transition from occ. beta-state:
            if statestr[transition[j][2]-1]=="2":
               sortstate="".join("u" if i==transition[j][2]-1 else c for i,c in enumerate(statestr))
            elif occ[1][transition[j][2]-1]==1:
               sortstate="".join("0" if i==transition[j][2]-1 else c for i,c in enumerate(statestr))
            else:
               print("  ",occ[0][transition[j][0]-1], transition[j][2])
         else:
            print(transition[j], statestr)
         if transition[j][1]!=0: #a transition to virt. alpha-state:
            if statestr[transition[j][1]-1]=="0":
               sortstate="".join("u" if i==transition[j][1]-1 else c for i,c in enumerate(sortstate))
            elif occ[1][transition[j][1]-1]==1:
               sortstate="".join("2" if i==transition[j][1]-1 else c for i,c in enumerate(sortstate))
         elif transition[j][3]!=0: #a transition to virt. beta-state:
            #print(j+1, transition[j], statestr[transition[j][3]-1], occ[1][transition[j][3]-1],occ[0][transition[j][3]-1])
            if statestr[transition[j][3]-1]=="0":
               sortstate="".join("d" if i==transition[j][3]-1 else c for i,c in enumerate(sortstate))
            elif occ[0][transition[j][3]-1]==1:
               sortstate="".join("2" if i==transition[j][3]-1 else c for i,c in enumerate(sortstate))

         output.write("%s"%sortstate)
         #print coefficient and weight
         sortstate=""
         output.write("     %16.10g         %16.10g\n"%( CIcoeff[i][j], CIcoeff[i][j]*CIcoeff[i][j]))
   output.close()

def printnorm(CIcoeff):
   nos=len(CIcoeff)
   trans=len(CIcoeff[0])
   for i in range(nos):
      norm=0
      for j in range(trans):
         norm+=CIcoeff[i][j]*CIcoeff[i][j]
      print(norm)

def printOrbitals(infile,outfile,state):
   output=open(outfile, 'a')
   files=open(infile, "r")
   inp=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
   files.close
   #for alpha-orbitals
   atemp=re.findall(
       b"(?<=DFT Final Alpha Molecular Orbital Analysis\n )[\w.=\+\- \n',^\"\d]+(?=DFT Final Beta)"
       , inp, re.M)[0]
   aMOvect=atemp.decode("utf-8").strip().split("Vector")
   anbf=len(aMOvect)-1 #because the first element is not an orbital vector
   aMOs, asort=getOrbitals(aMOvect[1] )
   #now, get the sorting and the first row to be printed
   #obtain the coefficients
   indMatrix=getCoefficients(aMOvect[1:])
   #print matrix in required format:
   if state=="f":
      output.write(u"\n\n final:")
   elif state=="i":
      output.write(u"\n\n initial:")
   else:
      output.write(u"\n\n unknown:")
   for n in range(0,anbf,5):
         output.write(u"\n\n alpha-Orbital")
         "                %i"
         for i in range(n,min(n+5,anbf)):
            output.write(u"                %i" %(i+1))
         output.write("\n\n")
         for j in range(len(aMOs)):
            output.write(u"  %2i    %s"%(j+1, aMOs[j]))
            for i in range(n,min(n+5,anbf)):
               output.write(u"  %0.10E " %(indMatrix[i][asort[j]])) #convert the repective element
            output.write("\n")
   output.close

   #for beta-orbitals
   btemp=re.findall(
       b"(?<=DFT Final Alpha Molecular Orbital Analysis\n )[\w.=\+\- \n',^\"\d]+(?=DFT Final Beta)"
       , inp, re.M)[0]
   bMOvect=btemp.decode("utf-8").strip().split("Vector")
   bnbf=len(bMOvect)-1 #because the first element is not an orbital vector
   bMOs, bsort=getOrbitals(bMOvect[1] )
   #test, if everything is consistent:
   if np.any(bsort!=asort):
      assert 1==2, "An error in the orbital-ordering occured. Alpha- and Beta-orbitals are sorted differently."
   assert bnbf==anbf, "Number of Basis functions differs between alpha and beta. This should never be the case."
   indMatrix=getCoefficients(bMOvect[1:])
   for n in range(0,bnbf,5):
         output.write(u"\n\n  beta-Orbital")
         "                %i"
         for i in range(n,min(n+5,bnbf)):
            output.write(u"                %i" %(i+1))
         output.write("\n\n")
         for j in range(len(bMOs)):
            output.write(u"  %2i    %s"%(j+1, bMOs[j]))
            for i in range(n,min(n+5,bnbf)):
               output.write(u"  %0.10E " %(indMatrix[i][bsort[j]])) #convert the repective element
            output.write("\n")
   output.close
   return asort #same as bsort

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
   
   roots=re.findall(b"(?<=Root )[ \d]+", inp)
   for i in range(len(roots)):
      roots[i]=int(roots[i])
   nos=np.max(roots)
   #cicoeff=re.findall(r"(?<=Dipole Oscillator Strength )[Occ. Virt. abE\-\+\d\n]+", inp)
   cicoeff=re.findall(b"(?<=Dipole Oscillator Strength )[Occ. Virt. abE\-\+\d\n alpha beta]+", inp)
   if len(cicoeff)>nos:
      cicoeff=cicoeff[-nos:] #this gives only last ci-vector
   elif len(cicoeff)<nos:
      assert 1==2, "an error occured. Not all roots found."
   for i in range(nos):
      CI=cicoeff[i].decode("utf-8").strip().split("\n")[2:] #split into lines and through first and last line away
      if i==0:
         if "Occ." in CI[-1]:
            noorb=len(CI)
         else:
            noorb=len(CI)-1
         #print(noorb)
         CIcoeff=np.zeros((nos,noorb))
         CItrans=np.zeros((noorb,4),dtype=int)
      for j in range(noorb):
         transition=CI[j].split()
         CIcoeff[i][j]=transition[-1]
         if transition[2]=="alpha":
            CItrans[j][0]=transition[1]
         elif transition[2]=="beta":
            CItrans[j][2]=transition[1]
         else:
            assert 1==2, "initial state unknown."
         if transition[7]=="alpha":
            CItrans[j][1]=transition[6]
         elif transition[7]=="beta":
            CItrans[j][3]=transition[6]
         else:
            assert 1==2, "final state unknown."
         #print(CItrans[j])
   noocc=int(max(np.max(CItrans[:].T[0]), np.max(CItrans[:].T[2])))
   nouno=int(max(np.max(CItrans[:].T[1]), np.max(CItrans[:].T[3])))-noocc
   return CIcoeff, CItrans, noocc, nouno, nos, noorb

def readOrbitals(infile):
   files=open(infile, "r")
   inp=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
   files.close
   atemp=re.findall(\
       b"(?<=DFT Final Alpha Molecular Orbital Analysis\n )[\w.=\+\- \n',^\"\d]+(?=DFT Final Beta)",
        inp, re.M)[-1]
   aMOvect=atemp.decode("utf-8").strip().split("Vector")
   anbf=len(aMOvect)-1 #because the first element is not an orbital vector
   aMOs, asort=getOrbitals(aMOvect[1] )
   #now, fet the sorting and the first row to be printed
   aoccupation=getOcc(aMOvect[1:])

   # repeat for beta-porbitals
   btemp=re.findall(b"(?<=DFT Final Beta Molecular Orbital Analysis\n )[\d\w .=\+\- \n',^\"]+(?=\n\n)", inp, re.M)[-1]
   #print(btemp)
   bMOvect=btemp.decode("utf-8").strip().split("Vector")
   bnbf=len(bMOvect)-1 
   bMOs, bsort=getOrbitals( bMOvect[1] )
   boccupation=getOcc(bMOvect[1:])

   #test, if asort==bsort
   if np.any(asort!=bsort):
      assert 1==2, "the sorting of alpha- and beta-orbitals doesn't coincide."
   occupation=[aoccupation, boccupation]
   MOs=[aMOs,bMOs]
   #print(MOs)
   return MOs, asort, occupation

def rOverlap(infile, sort):
   files=open(infile, "r")
   inp=mmap.mmap(files.fileno(), 0, prot=mmap.PROT_READ)
   files.close

   temp=re.findall(b"(?<=\n Begin overlap 1-e integrals\n ====================================\n)[\d.\w \+\- \n]+", inp)[-1]
   lines=temp.decode("utf-8").strip().split("\n")
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
   energy=float(re.findall(b"(?<=Total SCF energy =)[\d \-\.]+", inp)[-1])
   output.write("%15.10g\n"%energy)
   exen=re.findall(b"(?<=Root )[\w \d\.\-]+ ",inp)
   roots=[]
   for i in range(len(exen)):
      exen[i]=exen[i].decode("utf-8")
      if "." in exen[i]:
         roots.append(float(exen[i].strip().split()[0]))
   nos=int(np.max(roots))
   if "a.u." in exen[-1]:
      last=-1
   else:
      last=-2
   exen=exen[-nos+last:last] 
   for i in range(len(exen)):
      eorb=float(exen[i].strip().split()[-4])
      #assert eorb>0, "The state %s seems unconverged. Negative excitation energies occured." %(infile)
      output.write("%15.10g\n"%(eorb+energy))

def writePreamble(outfile, noocc,nouno, nbf, occi,occf, ftrans, itrans):
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
   output.write("FINALWF\n%d  %d\n\n"%(ftrans, ftrans))
   output.write("INITIALWF\n%d  %d\n\n"%(itrans, itrans))
   #output.write("INITIALWF\n1  1\n\n")
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

version=0.3
