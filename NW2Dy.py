#!/usr/bin/python3
import sys
from numpy import any
import unrest_rw2dy as rw


# changelog: 0.4:
# 1) use unrest_rw2dy instead rw2dy.
# 2) changed to python3
#
def main(argv=None):
   # evaluate input-arguments
   assert len(argv)==3, "three input files expected."
   fnfile=argv[0] # final file
   infile=argv[1] # initial file
   outfile=argv[2]
   
   # read first quantities from input-files
   fMO,fsort,occf=rw.readOrbitals(fnfile)
   MO,sort,occi=rw.readOrbitals(infile)
   #test, if everything is consistent so far
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
   #now, start writing to output-file
   rw.writePreamble(outfile, Inoocc,Inouno, dim, sum(occi[0]+occi[1]),sum(occf[0]+occf[1]), Ftrans, Itrans, Fnos, Inos)
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
   rw.printCI(outfile, FCIcoeff, FCitrans, sort,Fnoocc, Fnouno, Fnos,occf,Ftrans)
   output=open(outfile, 'a')
   output.write("INICI\nMULT=1\n")
   output.close()
   rw.printCI(outfile, ICIcoeff, ICitrans, sort,Inoocc, Inouno, Inos,occi,Itrans)

if __name__ == "__main__":
   main(sys.argv[1:])

version=0.4
