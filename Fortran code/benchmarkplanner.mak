# build the benchmark model of Khan and Thomas (2007)
# 03 - 05 April 2007

objects = benchmarkplanner.obj benchmarksim.obj erf.obj ppsplinefit3edit.obj

# complied using Intel Fortran Compiler Version 9.1  Build 20060927Z Package ID: W_FC_C_9.1.032
# using both the 64-bit (EMT64) and 32-bit versions
	  
compiler = ifort
optimise = O2
nobuild = c
stack = F512000000

benchmarkplanner.exe : $(objects) 
	$(compiler) /$(optimise) $(objects) /$(stack)
 
benchmarkplanner.obj : benchmarkplanner.f90 ppsplinefit3edit.obj kindset.obj
	$(compiler) benchmarkplanner.f90 /$(optimise) /$(nobuild)

benchmarksim.obj : benchmarksim.f90 kindset.obj
	$(compiler) benchmarksim.f90 /$(optimise) /$(nobuild)

include ppsplinefit3edit.mak

erf.obj : erf.f
	$(compiler) erf.f /$(optimise) /$(nobuild)

clean:
	del $(objects) kindset.mod ppsplinefit3edit.mod





