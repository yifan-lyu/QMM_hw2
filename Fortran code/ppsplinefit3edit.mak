# compile multivariate pp-spline module
# Aubhik Khan and Julia K. Thomas
# 03 April 2007

object = kindset.obj

# complied using Intel Fortran Compiler Version 9.1  Build 20060927Z Package ID: W_FC_C_9.1.032

compiler = ifort	
optimise = O2
nobuild = c

ppsplinefit3edit.obj : ppsplinefit3edit.f90 $(object)
	$(compiler) ppsplinefit3edit.f90 /$(optimise) /$(nobuild)

kindset.obj : kindset.f90
	$(compiler) kindset.f90 /$(optimise) /$(nobuild)





