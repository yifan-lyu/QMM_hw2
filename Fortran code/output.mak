# compile the output program for the inventory model of Khan and Thomas (2007)
# 07 July 2007

objects = outputinv.obj
compiler = ifort
optimise = Od
nobuild = c
stack = F512000000

# complied using Intel Fortran Compiler Version 10.0  Version 10.0    Build 20070426 Package ID: W_FC_P_10.0.025
# using both the 64-bit (EMT64) and 32-bit versions

outputinv.exe : $(objects) 
	$(compiler) /$(optimise) $(objects) /$(stack)
 
outputinv.obj : outputinv.f90 kindset.obj
	$(compiler) outputinv.f90 /$(optimise) /$(nobuild)

kindset.obj : kindset.f90
	$(compiler) kindset.f90 /$(optimise) /$(nobuild)



