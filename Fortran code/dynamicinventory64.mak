# build the (S,s) inventory model of Khan and Thomas (2007)
# 03 - 05 April 2007

objects = inventoriesedit.obj invworkshededit.obj invtoolshededit.obj \
	ppsplinefit3edit.obj betadistribution.obj kindset.obj betafunction.obj
objectsminor = ppsplinefit3edit.obj kindset.obj 	  
compiler = ifort
optimise = O2
nobuild = c
stack = F512000000

# include the appropriate files for ia32 or em64t compiler
# linkfiles = mkl_c.lib libguide.lib
# linkfiles = mkl_em64t.lib libguide.lib

linkfiles = mkl_c.lib libguide.lib

# use the path to the appropriate intel math kernel library for ia32 or em64t compiler
# libpath = libpath:"C:\Program Files\Intel\MKL\9.1.025\ia32\lib"
# libpath = libpath:"C:\Program Files\Intel\MKL\9.1.025\em64t\lib"

libpath = libpath:"C:\Program Files\Intel\MKL\9.1.025\ia32\lib"

# complied using Intel Fortran Compiler Version 10.0  Build 20070426 Package ID: W_FC_P_10.0.025
# using both the 64-bit (EMT64) and 32-bit versions

inventoriesedit.exe : $(objects) 
	$(compiler) /$(optimise) $(objects) /$(stack) /link $(linkfiles) /$(libpath)
 
inventoriesedit.obj : inventoriesedit.f90 invworkshededit.obj $(objectsminor)
	$(compiler) inventoriesedit.f90 /$(optimise) /$(nobuild)

invworkshededit.obj: invworkshededit.f90 invtoolshededit.obj $(objectsminor)
	$(compiler) invworkshededit.f90 /$(optimise) /$(nobuild)

invtoolshededit.obj : invtoolshededit.f90 $(objectsminor)
	$(compiler) invtoolshededit.f90 /$(optimise) /$(nobuild)

betadistribution.obj : betadistribution.f90 
	$(compiler) betadistribution.f90 /$(optimise) /$(nobuild)

include ppsplinefit3edit.mak

betafunction.obj : betafunction.for
	$(compiler) betafunction.for /$(optimise) /$(nobuild)






