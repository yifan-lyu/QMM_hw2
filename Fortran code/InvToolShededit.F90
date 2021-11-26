MODULE InvToolShededit

! Module of optimisation, regression and utility routines for the inventory model of Khan and Thomas (2007).

USE KINDSET
USE ppsplinefit3edit

IMPLICIT NONE

PRIVATE

PUBLIC:: InvM, InvSstar, CostFunction, InvPolicy, Setup, LinSpace, LogSpace

CONTAINS


!!!!!!!!!!!!!!!!!
!				!			
!	InvSstar	!
!				!
!!!!!!!!!!!!!!!!!

SUBROUTINE InvSstar(rknotsS, pq, precisionc, knotsS, csE1, Sbound, starscalar, Ea)


! Determine the target inventory stock using spline fitted to continuation value E1

USE KINDSET
use ppsplinefit3edit

IMPLICIT NONE

INTEGER(ik):: rknotsS, rknotsSint, s1
REAL(rk):: pq, precisionc, knotsS(rknotsS), Sbound(2), csE1(4,rknotsS-1), &
		   a, b, r, c, d, z, fc, fd, E1(1,1), price, star(1), Ea, &
		   starscalar 

INTENT(IN):: rknotsS, pq, precisionc, knotsS, csE1, Sbound
INTENT(OUT):: starscalar, Ea

rknotsSint = rknotsS - 2_ik

price = pq
s1 = 0_ik
z = 0.0_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							!	
!	Golden Section Search	!
!							!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Impose the constraint the m.LE.Sval hence Sf.GE.0 and apply golden section search
a = Sbound(1)
b = Sbound(2)

r = (3.0_rk - DSQRT(5.0_rk))/2.0_rk

! The c evaluation: Objective when m = c
c = a + r*(b-a) 

star(1) = c; CALL SPeval(csE1, knotsS, rknotsSint, 1, star, 1, 0, E1)
fc = -price*star(1) + E1(1,1); fc = -1.0_rk*fc

! The d evaluation
d = a + (1.0_rk - r)*(b-a)

star(1) = d; CALL SPeval(csE1, knotsS, rknotsSint, 1, star, 1, 0, E1)
fd = -price*star(1) + E1(1,1); fd = -1.0_rk*fd

DO
		  
	IF (DABS(d-c).LT.precisionc) THEN
		EXIT
	END IF
			                       
    s1 = s1 + 1_ik
    
	IF (fc.GE.fd) THEN
          
		z = c + (1.0_rk-r)*(b-c)
        ! [a c d b] <--- [c d z b]
        a = c
        c = d; fc = fd
		d = z 
	
		star(1) = d; CALL SPeval(csE1, knotsS, rknotsSint, 1, star, 1, 0, E1)
		fd = -price*star(1) + E1(1,1); fd = -1.0_rk*fd

	ELSE
		
		z = a + r*(d-a)
        ! [a c d b] <--- [a z c d]
        b = d
        d = c; fd = fc
        c = z; 
	
		star(1) = c; CALL SPeval(csE1, knotsS, rknotsSint, 1, star, 1, 0, E1)
		fc = -price*star(1) + E1(1,1); fc = -1.0_rk*fc

    END IF
    		           

END DO

starscalar = star(1)
Ea = -price*star(1) + E1(1,1)

END SUBROUTINE InvSstar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!													!!
!!		InvM: determination of M*, Sf=S1 - M*		!!
!!													!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Determine m given S1, production time inventory.

SUBROUTINE InvM(rknotsS, thetam, thetan, Atech, p, beta, precisiond, Apref, storage, storageh, knotsS, &
                csS, Sbound, zval, Sval, msol, nsol, cysol, E1sol)


USE KINDSET
use ppsplinefit3edit

IMPLICIT NONE

INTEGER(ik):: rknotsS, rknotsSint, s1
REAL(rk):: thetam, thetan, Atech, p, beta, precisiond, Apref, zval, &
		   Sval, knotsS(rknotsS), Sbound(2), csS(4,rknotsS-1), msol, nsol, g, &
		   a, b, r, c, d, z, fc, fd, V(1,1), Ve(2), m, Sf(1), E1sol, gz, cysol, rwsol, &
		   storage, storageh, fmsol, fnsol, fcysol, frwsol, fE1sol

INTENT(IN):: rknotsS, thetam, thetan, Atech, p, beta, precisiond, Apref, knotsS, &
             csS, Sbound, Sval, zval, storage, storageh
INTENT(OUT):: msol, nsol, E1sol, cysol

rknotsSint = rknotsS - 2_ik

gz = Atech*zval	

s1 = 0_ik

a=0.0_rk; b=0.0_rk; r=0.0_rk; c=0.0_rk; d=0.0_rk; z=0.0_rk; fc=0.0_rk; fd=0.0_rk
V=0.0_rk; Ve=0.0_rk; m=0.0_rk; Sf=0.0_rk; E1sol=0.0_rk; cysol=0.0_rk; rwsol=0.0_rk
msol=0.0_rk; nsol=0.0_rk; g=0.0_rk
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							!	
!	Golden Section Search	!
!							!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Impose the constraint the m.LE.Sval hence Sf.GE.0 and apply golden section search
a = Sbound(1)
b = Sbound(2)

r = (3.0_rk - DSQRT(5.0_rk))/2.0_rk

! The c evaluation: Objective when m = c
c = a + r*(b-a)

m = c; Sf(1) = Sval - m; CALL SPeval(csS, knotsS, rknotsSint, 1, Sf, 1, 0, V)
nsol = thetan*p*gz*(m**thetam)/Apref
nsol = nsol**(1.0_rk/(1.0_rk - thetan))
cysol = gz*(m**thetam)*(nsol**thetan) - storageh*storage*Sf(1)
rwsol = p*cysol - (1.0_rk - storageh)*p*storage*Sf(1) - Apref*nsol
fc = rwsol + beta*V(1,1); fc = -1.0_rk*fc

! The d evaluation
d = a + (1.0_rk - r)*(b-a)

m = d; Sf(1) = Sval - m; CALL SPeval(csS, knotsS, rknotsSint, 1, Sf, 1, 0, V)
nsol = thetan*p*gz*(m**thetam)/Apref
nsol = nsol**(1.0_rk/(1.0_rk - thetan))
cysol = gz*(m**thetam)*(nsol**thetan) - storageh*storage*Sf(1)
rwsol = p*cysol - (1.0_rk - storageh)*p*storage*Sf(1) - Apref*nsol
fd = rwsol + beta*V(1,1); fd = -1.0_rk*fd

DO
		  
	IF (DABS(d-c).LT.precisiond) THEN
		EXIT
	END IF
			                       
    s1 = s1 + 1_ik
               
    IF (fc.GE.fd) THEN
           
		z = c + (1.0_rk-r)*(b-c)
        ! [a c d b] <--- [c d z b]
        a = c
        c = d; fc = fd
		d = z 
		
		m = d; Sf(1) = Sval - m; CALL SPeval(csS, knotsS, rknotsSint, 1, Sf, 1, 0, V)
		nsol = thetan*p*gz*(m**thetam)/Apref
		nsol = nsol**(1.0_rk/(1.0_rk - thetan))
		cysol = gz*(m**thetam)*(nsol**thetan) - storageh*storage*Sf(1)
		rwsol = p*cysol - (1.0_rk - storageh)*p*storage*Sf(1) - Apref*nsol
		fd = rwsol + beta*V(1,1); fd = -1.0_rk*fd
                    
	ELSE
		
		z = a + r*(d-a)
        ! [a c d b] <--- [a z c d]
        b = d
        d = c; fd = fc
        c = z 
		
		m = c; Sf(1) = Sval - m; CALL SPeval(csS, knotsS, rknotsSint, 1, Sf, 1, 0, V)
		nsol = thetan*p*gz*(m**thetam)/Apref
		nsol = nsol**(1.0_rk/(1.0_rk - thetan))
		cysol = gz*(m**thetam)*(nsol**thetan) - storageh*storage*Sf(1)
		rwsol = p*cysol - (1.0_rk - storageh)*p*storage*Sf(1) - Apref*nsol 
		fc = rwsol + beta*V(1,1); fc = -1.0_rk*fc

    END IF
     
			           
END DO

msol = m; Sf(1) = Sval - msol; CALL SPeval(csS, knotsS, rknotsSint, 1, Sf, 1, 0, V)
nsol = thetan*p*gz*(m**thetam)/Apref
nsol = nsol**(1.0_rk/(1.0_rk - thetan))
cysol = gz*(m**thetam)*(nsol**thetan) - storageh*storage*Sf(1)
rwsol = p*cysol - (1.0_rk - storageh)*p*storage*Sf(1) - Apref*nsol
E1sol = rwsol + beta*V(1,1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							!
!	Exhaust S completely?	!
!							!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The possibility of fully exhausting the stock of intermediate goods.  Implies Sf = 0.0_rk hence no storage costs.

gz = Atech*zval

fmsol = Sval; Sf(1) = 0.0_rk; CALL SPeval(csS, knotsS, rknotsSint, 1, Sf, 1, 0, V)
fnsol = thetan*p*gz*(fmsol**thetam)/Apref
fnsol = fnsol**(1.0_rk/(1.0_rk - thetan))
fcysol = gz*(fmsol**thetam)*(fnsol**thetan)
frwsol = p*fcysol - Apref*fnsol
fE1sol = frwsol + beta*V(1,1)

IF (fE1sol.GE.E1sol) THEN
	msol = fmsol; nsol = fnsol; cysol = fcysol; E1sol = fE1sol
END IF

END SUBROUTINE InvM


!!	******************************************************************************************	!!
!!																								!!
!!																								!!
!!											CostFunction										!!
!!																								!!
!!																								!!
!!	******************************************************************************************	!!


FUNCTION CostFunction(zi, abeta, bbeta, masspoint, zibound, choice)

! 04.02.2002, generalized cost function
! phi = 0, A = 0 for our standard uniform distribution
! do not include choice variable if conditional expectation is needed

OPTIONAL:: choice
INTEGER(ik):: choice
REAL(rk):: zi, abeta, bbeta, zibound, masspoint, CostFunction, ziboundnew
INTENT(IN):: zi, abeta, bbeta, zibound, masspoint, choice

EXTERNAL BetaDistribution
 

IF (abeta.EQ.1.0_rk.AND.bbeta.EQ.1.0_rk) THEN

! uniform cost distribution generalized for a masspoint at an interior endpoint ziboundnew

ziboundnew = (1.0_rk - masspoint)*zibound

IF (PRESENT(choice)) THEN		! return G(zi) = probability z < zi
	IF (zi.LT.ziboundnew) THEN
		CostFunction = zi/zibound	
	ELSE
		CostFunction = 1.0_rk
	END IF
ELSE							! return conditional expectation	
	CostFunction = (zi**2.0_rk)/(2.0_rk*zibound)
	IF (zi.GE.ziboundnew) THEN
		CostFunction = CostFunction + masspoint*ziboundnew
	END IF
END IF

ELSE

! Generalized beta distribution

IF (PRESENT(choice)) THEN

	CALL BetaDistribution(abeta, bbeta, zi, zibound, masspoint, 0_ik, CostFunction)


ELSE

	CALL BetaDistribution(abeta, bbeta, zi, zibound, masspoint, 1_ik, CostFunction)

END IF

END IF

END FUNCTION CostFunction

!!	******************************************************************************************	!!
!!																								!!
!!																								!!
!!											InvPolicy											!!
!!																								!!
!!																								!!
!!	******************************************************************************************	!!


SUBROUTINE InvPolicy(T, Nx, Nz, ivec, Y, X, Beta, statistics)

! 19.12.2001
! 
! Subroutine to compute linear multivariate regression with an indicator 
! variable to sort separate equations.  There is a driver program.  
!
! Uses the Intel Math Kernel Library or equivalent library with BLAS

INTEGER(ik):: T, Nz, Nx, ivec(T), iz, places(Nz), info, it
REAL(rk):: Y(T), X(T,Nx), statistics(Nz,5), Beta(Nz, Nx + 1), XTX(Nx+1, Nx+1), &
		   XTY(Nx+1), P(Nx+1), ssr, ssd, R2, errormin, errormax, N	
REAL(rk), ALLOCATABLE:: Xs(:,:,:), Ys(:,:), estimates(:), residuals(:), &
					deviations(:) ! Xsi(:,:), Ysi(:)
INTENT(IN):: T, Nz, Nx, ivec, Y, X
INTENT(INOUT):: Beta, statistics

! sort the data
! each section (in the depth axis) of Xs and each column of Ys contains 
! the data for a distinct value of ivec arranged consecutively 
 
ALLOCATE(Ys(T, Nz), Xs(T, Nx+1, Nz))

Ys = 0.0_rk
Xs = 0.0_rk

places = 0_ik

DO iz = 1, Nz

	DO it = 1, T
		
		IF (ivec(it).EQ.iz) THEN
			places(iz) = places(iz) + 1_ik
			Xs(places(iz),1,iz) = 1.0_rk			! Add constant term			
			Xs(places(iz),2:Nx+1,iz) = X(it,1:Nx)	! Explanatory variables
			Ys(places(iz),iz) = Y(it)	 			! Dependent variable	
		END IF
	
	END DO

END DO

! regressions

DO iz = 1, Nz	! a separate regression for each value of ivec.
	
	XTX = MATMUL(TRANSPOSE(Xs(1:places(iz),1:Nx+1,iz)), Xs(1:places(iz),1:Nx+1,iz))
	XTY = MATMUL(TRANSPOSE(Xs(1:places(iz),1:Nx+1,iz)),Ys(1:places(iz),iz))

! Table 4-3, page 4-161, (340 of 1029 in mklman51.pdf)
! 
! CALL DGESV(Nx+1, nrhs, XTX, lda, ipiv, XTY, ldb, info)
! 
! INPUT 
! Nx+1 - order of XTX, rows of XTY
! nrhs - number of columns of XTY
! XTX, XTY the basic matrices of coefficients and constants, XTX*Beta = XTY
! lda - the first dimension of XTX, lda.GE.max(1,Nx+1)
! ldb - the first dimension of XTY, ldb.GE.max(1,Nx+1)
!
! OUTPUT
! XTX = PLU
! XTY = Beta, solution matrix
! P - permutation matrix
! info - output 0 if succesful, -i if i-th parameter was illegal, i if U(i,i) 
! if U(i,i) = 0 and the matrix U is singular so the solution could not be 
! completed after LU factorisation.

CALL DGESV(Nx+1, 1, XTX, Nx+1, P, XTY, Nx+1, info)

Beta(iz,1:Nx+1) = XTY

ALLOCATE(estimates(places(iz)), residuals(places(iz)), deviations(places(iz)))

N = DBLE(places(iz) - 2_ik)
estimates = MATMUL(Xs(1:places(iz),1:Nx+1,iz), Beta(iz,:))
residuals = Ys(1:places(iz),iz) - estimates
ssr = DOT_PRODUCT(residuals, residuals)
ssr = ssr/N
ssd = SUM(Ys(1:places(iz),iz))/DBLE(places(iz))
deviations =  Ys(1:places(iz),iz) - ssd
ssd = DOT_PRODUCT(deviations, deviations)
ssd = ssd/N
R2 = 1.0_rk - ssr/ssd

! Adjusted R2 
R2 = 1.0_rk - (1.0_rk - R2)*((N + 1.0_rk )/(N + 1.0_rk - DBLE(Nx)))

ssr = DSQRT(ssr)
errormin = MINVAL(ABS(residuals))
errormax = MAXVAL(ABS(residuals))

statistics(iz,1:5) = (/ssr, R2, errormin, errormax, N+2.0_rk/)	! contains ssr, R2 and error bounds

DEALLOCATE(estimates, residuals, deviations)

END DO

END SUBROUTINE InvPolicy


!!	******************************************************************************************	!!
!!																								!!
!!																								!!
!!					SetUp for Dynamic Program of the Inventory Model							!!
!!																								!!
!!																								!!
!!	******************************************************************************************	!!


SUBROUTINE SETUP(alpha, thetam, thetan, Atech, delta, abeta, bbeta, &
				 masspoint, zibound, outputcost, Apref, beta, storage, precisionv, &
				 precisionm, precisionS, precisionp, precisionq, NKM, KMbound, &
				 NS, NS0, Sbound, NSM, SMbound, pbounds, qbounds, &
				 simlength, Jmax, numbertry, dynamicfile, resultfile, shockfile, &
				 initialfile)  

INTEGER(ik):: NS, NS0, NKM, NSM, i, what, Jmax, simlength, numbertry
REAL(rk):: alpha, thetam, thetan, Atech, delta, beta, &
		   precisionv, precisionm, precisionS, precisionp, precisionq, Apref, &
		   Sbound(2), abeta, bbeta, masspoint, zibound, outputcost, &
		   pbounds(2), qbounds(2), KMbound(2), SMbound(2), switch, storage 
CHARACTER(30):: dynamicfile, resultfile, shockfile, initialfile

INTENT(OUT):: alpha, thetam, thetan, Atech, delta, abeta, bbeta, masspoint, &
			  zibound, Apref, beta, precisionv, precisionm, precisionS, precisionp, &
			  precisionq, NKM, KMbound, NS, NS0, Sbound, NSM, outputcost, &
			  SMbound, pbounds, qbounds, simlength, Jmax, &
			  dynamicfile, resultfile, shockfile, initialfile, storage

DO 

	WRITE(*,FMT = '(1X,A)', ADVANCE = 'NO') ' enter parameter file: '
	READ(*,*) dynamicfile
	
	
	OPEN (UNIT = 57, FILE = dynamicfile, STATUS = 'OLD', ACTION = 'READ', IOSTAT = what, ERR = 100)
	
	! The subroutine goes to line 100, right below, if there's a file opening error.  
	100	IF (what.NE.0_ik) THEN 
			WRITE(*,*) ' Hey, I think I could not locate ', dynamicfile
			WRITE(*,*) ' I got a file I/O error ', what
		ELSE
			EXIT
		END IF

END DO

READ(57,*)
READ(57,*)
READ(57,*)
READ(57,*) alpha										! capital's share in intermediate goods production
READ(57,*)
READ(57,*) thetam										! intermediate share in consumption production
READ(57,*)
READ(57,*) thetan										! labor's share in consumption production
READ(57,*) 
READ(57,*) Atech										! TFP term for cons. prod.
READ(57,*)
READ(57,*) delta										! depreciation	
READ(57,*)
READ(57,*) abeta										! beta distribution parameter a for G(zi)
READ(57,*)
READ(57,*) bbeta										! beta distribution parameter b for G(zi)
READ(57,*)
READ(57,*) zibound										! zi upper bar bound on adjustment costs
READ(57,*)
READ(57,*) masspoint									! masspoint on effective upper adjustment cost
READ(57,*)
READ(57,*) outputcost									! upper bound on uniformly distributed costs, G(zi)
READ(57,*)
READ(57,*) beta											! discount factor
READ(57,*)
READ(57,*) Apref										! preference term for leisure
READ(57,*)
READ(57,*) storage										! storage cost for inventories
READ(57,*)
READ(57,*) (Sbound(i), i = 1,2)							! knot endpoints on S
READ(57,*)
READ(57,*)  NS											! total number of knots for S
READ(57,*)
READ(57,*)  NS0											! total number of knots, including endpoints, on S	
READ(57,*)
READ(57,*) (KMbound(i), i = 1,2)						! knot endpoints on aggregate capital
READ(57,*)
READ(57,*)  NKM											! total number of knots, including endpoints, on KM	
READ(57,*)
READ(57,*) (SMbound(i), i = 1,2)						! knot endpoints on moment of S distribution
READ(57,*)
READ(57,*)  NSM											! total number of knots, including endpoints, on SM	
READ(57,*)
READ(57,*)  switch										! set to 1 for second order Taylor expansions
READ(57,*)  
READ(57,*)  precisionv									! precision term for value function convergence
READ(57,*)  
READ(57,*)  precisionm									! precision term for m choice
READ(57,*) 
READ(57,*)  precisionS									! precision term for S choice
READ(57,*) 
READ(57,*)  precisionp									! precision term for p choice
READ(57,*) 
READ(57,*)  precisionq									! precision term for q choice
READ(57,*)  
READ(57,*)  simlength									! simulation length
READ(57,*)
READ(57,*)  Jmax										! nonbinding maximum number of groups
READ(57,*)  
READ(57,*) (pbounds(i), i = 1,2)						! bisection endpoints on p
READ(57,*) 
READ(57,*) (qbounds(i), i = 1,2)						! gss bounds on q
READ(57,*) 
READ(57,*)  resultfile									! result file name
READ(57,*) 
READ(57,*)  shockfile									! exogenous shock file name
READ(57,*) 
READ(57,*)  initialfile									! initial steady state condition file
READ(57,*) 
READ(57,*)  numbertry									! initial steady state condition file

CLOSE (57)												! close datafile

END SUBROUTINE SETUP



!!	******************************************************************************************	!!
!!																								!!
!!																								!!
!!											LinSpace											!!
!!																								!!
!!																								!!
!!	******************************************************************************************	!!


SUBROUTINE LinSpace(lb, ub, gridnum, X)

INTEGER(ik):: gridnum, j1
REAL(rk):: lb, ub, X(gridnum), Y(gridnum-2_ik)

INTENT(IN):: lb, ub, gridnum
INTENT(OUT):: X

! This subroutine generates a vector which contains a grid 
! of n equally spaced elements between loweround and upper-
! bound.  Note this version produces a row vector.

DO j1 = 1_ik, gridnum - 2_ik
	Y(j1) = DBLE(j1)
END DO

! Note that I am multiplying an n-2x1 vector, Y, by a scalar 
! (which implies multiplying each element of the vector by 
! this scalar) and then adding another scalar (which implies 
! adding this scalar to each element of the vector).

Y = ((ub - lb)/DBLE(gridnum-1_ik))*Y + lb

X(1_ik) = lb
X(2_ik:gridnum-1_ik) = Y
X(gridnum) = ub

END SUBROUTINE LinSpace


!!	******************************************************************************************	!!
!!																								!!
!!																								!!
!!											LogSpace											!!
!!																								!!
!!																								!!
!!	******************************************************************************************	!!


SUBROUTINE LogSpace(lowerbound, upperbound, n, X)

USE KindSet

IMPLICIT NONE

INTEGER:: n, j1
REAL(rk):: lowerbound, upperbound, term0, lb, ub
REAL(rk):: X(n), Y(n-1)
INTENT(IN):: lowerbound, upperbound, n
INTENT(OUT):: X

term0 = DLOG(10.0_rk)
lb = DLOG(lowerbound)/term0
ub = DLOG(upperbound)/term0

DO j1 = 1, n - 2
	Y(j1) = DBLE(j1)
END DO

Y = ((ub - lb)/dble(n-1))*Y + lb

X(1) = lb
X(2:n-1) = Y
X(n) = ub

X = 10.0_rk**X

END SUBROUTINE LogSpace

END MODULE InvToolShededit