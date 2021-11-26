# compile the steady state of the inventory model of Khan and Thomas (2007)
# 05 April 2007

objects = steadystate.obj betadistribution.obj betafunction.obj ppsplinefit3edit.obj
compiler = ifort
optimise = O2
nobuild = c

# complied using Intel Fortran Compiler Version 10.0  Version 10.0    Build 20070426 Package ID: W_FC_P_10.0.025
# using both the 64-bit (EMT64) and 32-bit versions

steadystate.exe : $(objects) 
	$(compiler) /$(optimise) $(objects) 
 
steadystate.obj : steadystate.f90 betadistribution.obj kindset.obj ppsplinefit3edit.obj
	$(compiler) steadystate.f90 /$(optimise) /$(nobuild)

betadistribution.obj : betadistribution.f90 betafunction.obj kindset.obj
	$(compiler) betadistribution.f90 /$(optimise) /$(nobuild)

include ppsplinefit3edit.mak

betafunction.obj : betafunction.for
	$(compiler) betafunction.for /$(optimise) /$(nobuild)






