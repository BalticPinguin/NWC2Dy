#!/usr/bin/python3.4
import sys
from numpy import any
import numpy as np
import unrest_rw2dy as rw


# changelog: 0.4:
# 1) use unrest_rw2dy instead rw2dy.
# 2) changed to python3
# 3) introduced cutoff-energies -> simple implementation
#
def main(argv=None):
   # evaluate input-arguments
   assert len(argv)==3, "three input files expected."
   fnfile=argv[0] # final file
   infile=argv[1] # initial file
   outfile=argv[2]
   
   # read first quantities from input-files
   fMO,fsort,occf,enf=rw.readOrbitals(fnfile)
   MO,sort, occi, eni=rw.readOrbitals(infile)
 
   lowen=5000000
   highen=500000
   HOMO=np.max(np.where(occf[0]==1)) #HOMO==FNOOCC !?
   try:
      low=np.max(np.where(enf[0][HOMO]-enf[0]>lowen))-1
   except ValueError: # lowen larger thean all orbital energies:
      low=0
   try:
      high=np.min(np.where(enf[0]-enf[0][HOMO]>highen))
   except ValueError: # lowen larger thean all orbital energies:
      high=len(occf[0])

#   print(lowboarder)
#   print(highboarder)
#
#   print("final state\n alpha-occ alpha-en   beta-occ beta-en")
#   for i in range(len(occf[0])):
#      print(occf[0][i], enf[0][i], " ", occf[1][i], enf[1][i])
#   print("\n\ninitial state\n alpha-occ alpha-en   beta-occ beta-en")
#   for i in range(len(occf[0])):
#      print(occi[0][i], eni[0][i], " ", occi[1][i], eni[1][i])
   
   #frozen=low+len(occf[0])-high
   frozen=0 #since low=...-1 to make it index-corrected.
   if low==0:
      if high>0:
         occf[0]=occf[0][:high] #truncate occupation vectors accoringly.
         occf[1]=occf[1][:high]
         occi[0]=occi[0][:high]
         occi[1]=occi[1][:high]
   #else:  keep it as it is!
   else:
      if high>0:
         occf[0]=occf[0][low:high] #truncate occupation vectors accoringly.
         occf[1]=occf[1][low:high]
         occi[0]=occi[0][low:high]
         occi[1]=occi[1][low:high]
      else:
         occf[0]=occf[0][low:]
         occf[1]=occf[1][low:]
         occi[0]=occi[0][low:]
         occi[1]=occi[1][low:]   #test, if everything is consistent so far
   if MO!=fMO:
      error=open("error.out", 'w')
      for i in range(len(MO)):
         error.write("%.8g   %.8g\n"%(MO[i], fMO[i]))
      error.close()
      assert 1==0, "the MO-vectors don't coinciede. See file error.out"
   if any(sort!=fsort):
      error=open("error.out", 'w')
      for i in range(len(sort)):
         error.write("%.8g   %.8g    %8g\n"%(sort[i], fsort[i], sort[i]-fsort[i]))
      error.close()
      assert 1==0, "the The sorting-vectors don't coinciede. See file error.out"
   ov, dim=rw.rOverlap(infile, sort)
   assert sum(occi[0]+occi[1])==sum(occf[0]+occf[1])+1,\
      "the numbers of initial/final electrons seem to be wrong:"+\
      " initial: %i final: %i" %(sum(occi[0]+occi[0]), sum(occf[0]+occf[1]))

 ######## this is the alternative way: read CI-vectors from extra files
 #  CIfile=getFile("CI-coefficients")
 #  CIcoeff, Citrans, noocc, nouno, nos=readCI(CIfile)
 ########

   #read CI-vectors
   FCIcoeff, FCitrans, Fnoocc, Fnouno, Fnos,Ftrans=rw.readCI2(fnfile)
   ICIcoeff, ICitrans, Inoocc, Inouno, Inos,Itrans=rw.readCI2(infile)
   #print(frozen, high, low, Fnoocc, Fnouno)

   #renormalize vectors (due to cutoff)
   FCIcoeff=rw.normalise(FCIcoeff)
   ICIcoeff=rw.normalise(ICIcoeff)

   #now, start writing to output-file
   rw.writePreamble(outfile, Inoocc,Inouno, dim, frozen, sum(occi[0]+occi[1]), 
                            sum(occf[0]+occf[1]), Ftrans, Itrans, Fnos, Inos)
   rw.printOrbitals(fnfile,outfile,2,"f")
   rw.printOrbitals(infile,outfile,1,"i")
   rw.wOverlap(ov,dim,sort,outfile)
   #printnorm(CIcoeff) # this is for debugging only
   output=open(outfile, 'a')
   output.write("FINEN \n")
   rw.rwenergy(output,fnfile)
   output.write("\nINIEN \n")
   rw.rwenergy(output,infile)
   output.write("\n\nFINCI\nMULT=2\n")
   output.close()
   rw.printCI(outfile, FCIcoeff, FCitrans, sort,Fnoocc, Fnouno, Fnos, occf, Ftrans)
   output=open(outfile, 'a')
   output.write("INICI\nMULT=1\n")
   output.close()
   rw.printCI(outfile, ICIcoeff, ICitrans, sort,Inoocc, Inouno, Inos,occi,Itrans)

if __name__ == "__main__":
   main(sys.argv[1:])

version=0.4
