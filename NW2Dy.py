#!/usr/bin/python
import sys
import unrest_rw2dy as rw


# changelog: 0.4:
# 1) use unrest_rw2dy instead rw2dy.
#
def main(argv=None):
   assert len(argv)==3, "three input files expected."
   fnfile=argv[0] # final file
   infile=argv[1] # initial file
   outfile=argv[2]
   fMO,fsort,occf=rw.readOrbitals(fnfile)
   MO,sort,occi=rw.readOrbitals(infile)
   #test, if everything is fine:
   if MO!=fMO:
      error=open("error.out", 'w')
      for i in range(len(MO)):
         error.write("%.8g   %.8g"%(MO[i], fMO[i]))
      error.close()
      assert 1==0, "the MO-vectors don't coinciede. See file error.out"
   if any(sort!=fsort):
      error=open("error.out", 'w')
      for i in range(len(MO)):
         error.write("%.8g   %.8g"%(sort[i], fsort[i]))
      error.close()
      assert 1==0, "the The sorting-vectors don't coinciede. See file error.out"
   ov, dim=rw.rOverlap(infile, sort)
   if sum(occi[0]+occi[1])!=sum(occf[0]+occf[1])+1:
      assert 2==1, "the numbers of initial/final electrons seem to be wrong:"\
      " initial: %i final: %i" %(sum(occi[0]+occi[0]), sum(occf[0]+occf[1]))
 
 ######## this is the alternative way: read CI-vectors from extra files
 #  CIfile=getFile("CI-coefficients")
 #  CIcoeff, Citrans, noocc, nouno, nos=readCI(CIfile)
 ########

   FCIcoeff, FCitrans, Fnoocc, Fnouno, Fnos=rw.readCI2(fnfile)
   ICIcoeff, ICitrans, Inoocc, Inouno, Inos=rw.readCI2(infile)

   #now, start writing to output-file
   #writePreamble(outfile, Inoocc,Inouno, dim)
   rw.writePreamble2(outfile, Inoocc,Inouno, dim, sum(occi[0]+occi[1]),sum(occf[0]+occf[1]))
   rw.printOrbitals(infile,outfile)
   rw.wOverlap(ov,dim,sort,outfile)
   #printnorm(CIcoeff)
   output=open(outfile, 'a')
   output.write("\nFINEN \n")
   rw.rwenergy(output,fnfile)
   output.write("\nINIEN \n")
   rw.rwenergy(output,infile)
   output.close()
   rw.printCI(outfile, FCIcoeff, FCitrans, sort,Fnoocc, Fnouno, Fnos,occf)
   rw.printCI(outfile, ICIcoeff, ICitrans, sort,Inoocc, Inouno, Inos,occi)

if __name__ == "__main__":
   main(sys.argv[1:])

version=0.4
