#!/usr/bin/python
import sys
import rw2dy as rw

def main(argv=None):
   assert len(argv)==3, "three input files expected."
   fnfile=argv[0]
   infile=argv[1]
   outfile=argv[2]
   MO,sort,occi=rw.readOrbitals(infile)
   fMO,fsort,occf=rw.readOrbitals(fnfile)
   if MO!=fMO:
      print "foo"
   if any(sort!=fsort):
      print "bar"
   ov, dim=rw.rOverlap(infile, sort)
   #CIfile=getFile("CI-coefficients")
   #CIcoeff, Citrans, noocc, nouno, nos=readCI(CIfile)
   FCIcoeff, FCitrans, Fnoocc, Fnouno, Fnos=rw.readCI2(fnfile)
   ICIcoeff, ICitrans, Inoocc, Inouno, Inos=rw.readCI2(infile)
   #test, if everything is fine:

   #now, start writing to output-file
   #writePreamble(outfile, Inoocc,Inouno, dim)
   rw.writePreamble2(outfile, Inoocc,Inouno, dim, sum(occi),sum(occf))
   rw.printOrbitals(infile,outfile)
   rw.wOverlap(ov,dim,sort,outfile)
   #printnorm(CIcoeff)
   output=open(outfile, 'a')
   output.write("\nFINEN \n")
   rw.rwenergy(output,fnfile)
   output.write("\nINIEN \n")
   rw.rwenergy(output,infile)
   output.close()
   #print occf
   #print occi
   rw.printCI(outfile, FCIcoeff, FCitrans, sort,Fnoocc, Fnouno, Fnos,occf)
   rw.printCI(outfile, ICIcoeff, ICitrans, sort,Inoocc, Inouno, Inos,occi)

if __name__ == "__main__":
   main(sys.argv[1:])

version=0.3
