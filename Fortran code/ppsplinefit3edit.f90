MODULE PPSplineFit3edit

! This module contains subroutines that implement piecewise polynomial cubic spline 
! interpolation.  Interpolation of functions of arbitrary dimension is possible, and 
! the algorithm computes function values as well as first and second derivatives.  
!   
! For univariate interpolation, complete, natural and not-a-knot endpoint conditions 
! are allowed.  Multivariate interpolation is implemented using the algorithm described in 
!
! Johnson, S. A. (1989) Spline approximation in discrete dynamic programming with application
! to stochastic multi-reservoir systems Unpublished dissertation (Cornell, Ithaca, NY).
!
! and requires the not-a-knot endpoint  condition.  
! 
! Aubhik Khan and Julia K. Thomas, 04.09.2003, 23.11.2003.
!
! Edited for distribution on 18.02.2007

USE KindSet

IMPLICIT NONE

PRIVATE

PUBLIC SPLHS, SPpp, SPEval, SPFitA0, SPFitA, SPFitB, FastSplinEval

CONTAINS


! ************************************************************************************ !
 
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!!					!!
!!		SPLHS		!!	
!!					!!
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!

! This subroutine within a module contains, within itself, a subroutine

SUBROUTINE SPLHS(knots,r,indicator,L,U,dtau)

! Dimensions
!
!				r		number of interior knots
!				knots	1 x (r + 2)
!				L		1 x (r-1)
!				U		2 x r
!				dtau	r+1 x 1

INTEGER, PARAMETER:: z = 1, pind = 2, lind = 3, uind = 1
INTEGER(ik), INTENT(IN):: r
INTEGER(ik):: j
REAL(rk), INTENT(IN):: knots(z+r+1,1)
CHARACTER(20), INTENT(IN):: indicator
REAL(rk), INTENT(OUT):: L(1,r-1), U(2,r), dtau(r+1,1)
REAL(rk)::  T(3,r), ratio0, ratior, NAK11, NAK12, NAKrr, NAKrrm1

! Note: must specifiy dimension of dtau below.  
dtau(1:r+1,1) = knots(z+1:z+r+1,1) - knots(z:z+r,1)

T = 0.0_rk

IF (indicator(1:5).EQ.'not-a') THEN			! Not-a-knot endpoint condition

   ratio0	= dtau(z,1)/dtau(z+1,1)
   ratio0	= ratio0*dtau(z,1)
   NAK11	= -dtau(z+1,1) + ratio0
   NAK12	= ratio0   
   ratior	= dtau(z+r,1)/dtau(z+r-1,1)
   ratior	= ratior*dtau(z+r,1)
   NAKrr	= -dtau(z+r-1,1) + ratior
   NAKrrm1	= ratior;

ELSEIF (indicator(1:5).EQ.'natur') THEN		! Natural endpoint condition

	NAK11 = -0.5_rk*dtau(z+1,1)
	NAK12 = 0.0_rk
	NAKrr = -0.5_rk*dtau(z+r-1,1)
	NAKrrm1 = 0.0_rk
ELSE

   NAK11	= 0.0_rk
   NAK12	= 0.0_rk
   NAKrr	= 0.0_rk
   NAKrrm1	= 0.0_rk
END IF
   
! This matrix is entered differently from T in the notes, for 
! computational reasons.  In any case, note that lind is the lower index
! equivalent to T(j-1,j), pind is the anchor T(j,j) and uind is the 
! upper index T(j+1,j).  Endpoint conditions for natural or not-a-knot 
! splines appear in NAK11, NAK12, NAKrr and NAKrrm1.    

T(pind,1) = NAK11 + 2.0_rk*(dtau(z,1) + dtau(z+1,1))
T(uind,1) = NAK12 + dtau(z,1)

T(lind,2) = dtau(z+2,1)
T(uind,2) = dtau(z+1,1)
T(pind,2) = 2.0_rk*(T(lind,2) + T(uind,2))

DO j = 3,r-1,1
   T(lind,j) = dtau(z+j,1)
   T(uind,j) = T(lind,j-1)
   T(pind,j) = 2.0_rk*(T(lind,j) + T(uind,j))
END DO

T(pind,r) = NAKrr + 2.0_rk*(dtau(z+r,1) + dtau(z+r-1,1))
T(lind,r) = NAKrrm1 + dtau(z+r,1)

CALL TridiagLUfactor(T, r, pind, uind, lind, L, U)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!
!						!	
!	TriDiagLUfactor		!
!						!
!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TridiagLUfactor(T, n, pind, uind, lind, L, U)

INTEGER, INTENT(IN):: n, pind, uind, lind
INTEGER(ik):: i
REAL(rk), INTENT(IN):: T(3,r)
REAL(rk), INTENT(OUT):: L(1,r-1), U(2,r)

L = 0.0_rk
U = 0.0_rk
U(uind,1:n-1) = T(uind,1:n-1)
U(pind,1) = T(pind,1)

DO i = 2,n-1,1
   L(1,i-1) = T(lind,i)/U(pind,i-1)
   U(pind,i) = T(pind,i) - L(1,i-1)*U(uind,i-1)
END DO

L(1,n-1) = T(lind,n)/U(pind,n-1)
U(pind,n) = T(pind,n) - L(1,n-1)*U(uind,n-1)

END SUBROUTINE TridiagLUfactor

END SUBROUTINE SPLHS

! ************************************************************************************ !

!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!!					!!
!!		SPpp		!!
!!					!!
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!

! Dimensions
!
!				m		number of univariate functions being splined
!				fvals	(r+2) x m
!				SVEC	2 x m
!				c		4 x (r+1)*m

SUBROUTINE SPpp(fvals, r, m, L, U, dtau, indicator, SVEC, c)

INTEGER(ik), PARAMETER:: z = 1_ik
INTEGER(ik), INTENT(IN):: r, m
REAL(rk), INTENT(IN):: fvals(z+r+1,m), SVEC(2,m), L(1,r-1), U(2,r), dtau(r+1,1)
REAL(rk):: B(r,m), y(r,m), s(r,m), df(z+r,m)
REAL(rk), INTENT(OUT):: c(4,(r+z)*m)
CHARACTER(20), INTENT(IN):: indicator

! Program accepts function values at a grid of r+2 knot points for m functions,
! the LU factorisation of the LHS of the system of equations which determine
! the slope coefficients and are invariant to function values, and the similarly 
! invarant function of the difference of the knots, dtau.  It also requires a 
! choice of 'not-a-knot' or complete spline endpoint conditions in indicator and, 
! if complete is chosen, the endpoint slopes for each of the m splines being 
! simultaneously determined.  

! Determine B, the RHS to the m system of equations and df, divided difference of 
! function values.
CALL SPRHS(dtau, r, m, fvals, indicator, SVEC, B, df)

! Determine y from B as the first step in solving the system of equations
CALL SPLsolve(L,B,r,m,y)

! Determine r interior slopes s(i,j), i = 1,...,r, j = 1,...,m.  s(:,j) describes 
! the r interior slopes for a given spline of which there are m.
CALL SPUsolve(U,y,r,m,s)

! Use s, fvals, dtau and df to retrieve coefficients for the m splines in c.
CALL SPppcoeff(fvals, dtau, r, m, s, df, indicator, SVEC, c)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!
!					!	
!	SPRHS (SPpp)	!
!					!
!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SPRHS(dtau, r, m, fvals, indicator, SVEC, f, df)

INTEGER(ik), PARAMETER:: z = 1
INTEGER(ik), INTENT(IN):: r, m
REAL(rk), INTENT(IN):: dtau(r+z,m), fvals(z+r+1,m), SVEC(2,m)
CHARACTER(20), INTENT(IN):: indicator
REAL(rk), INTENT(OUT):: f(r,m), df(z+r,m)
INTEGER(ik):: j
REAL(rk):: ratio0(m), NAK1(m), ratior(m), NAKr(m)

! This program is able to simultaneously set up the RHS of the 
! slope equations which determine the slope terms for each of m 
! splines defined upon the same r+2 set of knots.  This program 
! too was extended to accomodate natural splines on 20.03.2002.
! Additional columns of f were introduced to allow for 
! multiple simultaneous spline fits on the same set of knots.

df = fvals(z+1:z + r + 1,:) - fvals(z:z + r,:)
df = df/dtau

f = 0.0_rk

IF (indicator(1:5).EQ.'not-a') THEN							
   
   ! A not-a-knot spline
   ratio0 = dtau(z,:)/dtau(z+1,:)
   ratio0 = ratio0*dtau(z,:)
   NAK1 = -2.0_rk*dtau(z+1,:)*df(z,:) + 2.0_rk*ratio0*df(z+1,:)
   
   ratior = dtau(r+z,:)/dtau(r+z-1,:)
   ratior = ratior*dtau(r+z,:)
   NAKr = 2.0_rk*ratior*df(r+z-1,:) - 2.0_rk*dtau(r+z-1,:)*df(r+z,:)
   
ELSEIF(indicator(1:5).EQ.'natur') THEN

	! For a natural spline it is assumed that SVEC(1) = mu0 and SVEC(2) = mur,
	! the second derivative values at the endpoints.
	NAK1 = -1.5_rk*dtau(z+1,:)*df(z,:) + dtau(z,:)*dtau(z+1,:)*SVEC(1,:)/4.0_rk
	NAKr = -1.5_rk*dtau(z+r-1,:)*df(z+r,:) - dtau(z+r-1,:)*dtau(z+r,:)*SVEC(2,:)/4.0_rk

ELSE
      
	! For a complete spline it is assumed that SVEC(1) = s0 and SVEC(2) = sr+1,
	! the slopes at the endpoints.
    NAK1 = -dtau(z+1,:)*SVEC(1,:)
    NAKr = -dtau(z+r-1,:)*SVEC(2,:)
   
END IF

f(1,:) = dtau(z+1,:)*df(z,:) + dtau(z,:)*df(z+1,:)
f(1,:) = 3.0_rk*f(1,:) + NAK1
f(r,:) = dtau(r+z,:)*df(r+z-1,:) + dtau(r+z-1,:)*df(r+z,:)
f(r,:) = 3.0_rk*f(r,:) + NAKr

DO j = 2,r-1,1
   f(j,:) = dtau(z+j,:)*df(z + j - 1,:) + dtau(z + j - 1,:)*df(z + j,:);
   f(j,:) = 3.0_rk*f(j,:)
END DO

END SUBROUTINE SPRHS


!!!!!!!!!!!!!!!!!!!!!!!!!
!						!	
!	 SPLsolve (SPpp)	!
!						!
!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SPLsolve(L,B,r, m, y)

INTEGER(ik), INTENT(IN):: r, m 
REAL(rk), INTENT(IN):: L(1,r-1), B(r,m)
REAL(rk), INTENT(OUT):: y(r,m)
INTEGER(ik):: i

y = 0.0_rk
y(1,:) = B(1,:);

DO i = 2, r, 1
   y(i,:) = B(i,:) - L(1,i-1)*y(i-1,:);
END DO


END SUBROUTINE SPLsolve


!!!!!!!!!!!!!!!!!!!!!
!					!
!	SPUsolve (SPpp)	!
!					!
!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE SPUsolve(U,y,r,m,x)

INTEGER(ik), PARAMETER:: pind = 2, uind = 1
INTEGER(ik), INTENT(IN):: r, m
INTEGER(ik):: i
REAL(rk), INTENT(IN):: U(2,r), y(r,m)
REAL(rk), INTENT(OUT):: x(r,m)

x =  0.0_rk;
x(r,:) = y(r,:)/U(pind,r);

DO i = r-1,1,-1
   x(i,:) = (y(i,:) - U(uind,i)*x(i+1,:))/U(pind,i);
END DO

END SUBROUTINE SPUsolve


!!!!!!!!!!!!!!!!!!!!!!!!!
!						!
!	SPppcoeff (SPpp)	!
!						!
!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE SPppcoeff(fvals, dtau, r, m, s, df, indicator, SVEC,c)

INTEGER(ik), PARAMETER:: z = 1
INTEGER(ik), INTENT(IN):: r, m
INTEGER(ik):: nc
REAL(rk), INTENT(IN):: fvals(z+r+1,m), dtau(r+z,m) , s(r,m), df(z+r,m), SVEC(2,m)
REAL(rk), INTENT(OUT):: c(4,(r+z)*m)
REAL(rk):: dtau2(r+z,m), dratios0(m), dratiosrp1(m), s0(m), srp1(m), sadj(z+r+1,m), &
	junk(1,(r+z)*m), junk2(1,(r+z)*m)
CHARACTER(20), INTENT(IN):: indicator

! This program is able to fit m univariate splines simultaneously, each defined 
! on the same r+2 set of knots.  The splines are stored consecutively in c, the j-th
! spline has r+1 polynomial pieces described by the columns of coefficients in from 
! (j-1)*(r+1) + 1 to j*(r+1) in c.

! Part 1) Use endpoint conditions to retrieve endpoint slopes

dtau2 = dtau*dtau
nc = (r+z)*m										! for each of m splines, r+1 polynomial pieces,  
c = 0.0_rk											! each a column of 4 coefficients within c	

IF (indicator(1:5).EQ.'not-a' ) THEN				! not-a-knot endpoint condition
   
   dratios0 = dtau(z,:)/dtau(z+1,:)					
   dratios0 = dratios0*dratios0
   dratiosrp1 = dtau(z+r,:)/dtau(z+r-1,:)
   dratiosrp1 = dratiosrp1*dratiosrp1
   
   s0 = (-1.0_rk + dratios0)*s(1,:)
   s0 = s0 + dratios0*s(2,:) + 2.0_rk*df(z,:) - 2.0_rk*dratios0*df(z+1,:)
   
   srp1 = (-1.0_rk + dratiosrp1)*s(r,:)
   srp1 = srp1 + dratiosrp1*s(r-1,:) + 2.0_rk*df(z+r,:) - 2.0_rk*dratiosrp1*df(z+r-1,:)

ELSEIF (indicator(1:5).EQ.'natur' ) THEN
	
	s0 = 1.5_rk*df(z,:) - 0.5_rk*s(1,:) - dtau(z,:)*SVEC(1,:)/4.0_rk
	srp1 = 1.5_rk*df(z+r,:) - 0.5_rk*s(r,:) + dtau(z+r,:)*SVEC(2,:)/4.0_rk
   
ELSE												! complete spline 		
   
   s0 = SVEC(1,:)
   srp1 = SVEC(2,:)
   
END IF

sadj(z,1:m) = s0
sadj(z+1:r+1,1:m) = s
sadj(z+r+1,1:m) = srp1

! Part 2) retrieve coefficients for each polynomial

junk = RESHAPE((fvals(z:r+z,:)),(/1,nc/))			
c(1,:) = junk(1,:)

junk  = RESHAPE(sadj(z:r+z,:),(/1,nc/))				
c(2,:) = junk(1,:)

junk = RESHAPE((df(z:r+z,:)/dtau(z:r+z,:)),(/1,nc/))
junk2 = RESHAPE((sadj(z:r+z,:)/dtau(z:r+z,:)),(/1,nc/))
junk = 3.0_rk*junk - 2.0_rk*junk2
junk2 = RESHAPE((sadj(z+1:r+z+1,:)/dtau(z:r+z,:)),(/1,nc/))
junk = junk - junk2
c(3,:) = junk(1,:)										

junk = -1.0_rk*RESHAPE((df(z:r+z,:)/dtau2(z:r+z,:)),(/1,nc/))
junk2 = RESHAPE((sadj(z:r+z,:)/dtau2(z:r+z,:)),(/1,nc/))
junk = 2.0_rk*junk + junk2
junk2 = RESHAPE((sadj(z+1:r+z+1,:)/dtau2(z:r+z,:)),(/1,nc/))
junk = junk + junk2
c(4,:) = junk(1,:)


END SUBROUTINE SPppcoeff

END SUBROUTINE SPpp


! ************************************************************************************ !

!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!!					!!
!!  	SPeval		!!
!!					!!
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!

! Dimensions
!
!				zoom	number of points being evaluated
!				order	0,1 or 2 for function value, first or second derivative
!				x		vector of length zoom
!				h		(order+1) x (zoom*m)

SUBROUTINE SPeval(c, knots, r, m, x, zoom, order, h)

INTEGER(ik):: r, zoom, order, knotmin, knotmax, knothigh, knotlow, i, knotloc(zoom), m, &
	knotplace, j
REAL(rk):: x(zoom), c(4, (r+1)*m), knots(r+2), y(zoom), h(order + 1, zoom*m)

INTENT(IN):: c, knots, r, m, x, zoom, order
INTENT(OUT):: h

! There are m separate splines stored column consecutively in c, which is 4 x (r+1)*m
! hence r is the number of interior knots and r+1 the number of polynomial pieces for
! each spline.  x is a row vector of 1 x zoom points at which each spline is to be 
! evaluated.  

! Part 1) locate the appropriate knot

! We use a bisection search to locate, on the integer grid spaces of polynomial spline 
! pieces, the location of tau(i) such that x(i) - tau(i) is the smallest positive 
! number.  

knotmin = 1_ik;			! left endpoint knot plus one
knotmax = r+2_ik;		! right endpoint knot
knothigh = knotmax;		! current upper bound on corresponding knot place
knotlow = knotmin;

DO i = 1, zoom											! must locate polynomial piece for each value of x

    knothigh = knotmax											
	knotlow = knotmin
		
  DO 

  IF (knothigh - knotlow.LE.1_ik) THEN					! exit indefinite do loop (while loop)
	EXIT
  END IF      

     knotplace = (knothigh + knotlow)/2_ik				! bisect

     IF (knots(knotplace).GE.x(i)) THEN					! is this tau(i) >= x(i)?
        knothigh = knotplace;							! yes, hence it is the upper bound on the bracket
     ELSE
        knotlow = knotplace;							! no, hence it is the lower bound on the bracket
     END IF
	      
  END DO													
     
  knotloc(i) = knotlow;									! the polynomial piece for x(i)

END DO

! Part 2) return spline value using corresponding polynomial piece

DO i = 1,zoom											! for each value of x
	y(i) = x(i) - knots(knotloc(i))
END DO


DO j = 1,m												! for each univariate spline

	IF (order.EQ.0) THEN								! return function values only

		h(1, (j-1)*zoom + 1: j*zoom) = c(4,knotloc)
 
		DO i = 3, 1, -1
			h(1,(j-1)*zoom + 1: j*zoom) = h(1,(j-1)*zoom + 1: j*zoom)*y + c(i,knotloc)
		END DO

	ELSE												! return function and first derivative values

		h(1, (j-1)*zoom + 1: j*zoom) = c(4,knotloc)		! function values in first row of h

		h(2, (j-1)*zoom + 1: j*zoom)  = 0.0_rk			! first derivative values in second

		DO i = 3, 1, -1
			h(1,(j-1)*zoom + 1: j*zoom) = h(1,(j-1)*zoom + 1: j*zoom)*y + c(i, knotloc)
			h(2,(j-1)*zoom + 1: j*zoom) = h(2,(j-1)*zoom + 1: j*zoom)*y + DBLE(i)*c(i+1, knotloc)
		END DO
   
		IF (order.GT.1) THEN							! add second derivatives in third row
			h(3,(j-1)*zoom + 1: j*zoom) = 2.0_rk*c(3,knotloc) + 6.0_rk*c(4,knotloc)*y
		END IF
      
	END IF

	knotloc = knotloc + (r+1)							! shift up to next set of spline coefficients	

END DO 

END SUBROUTINE SPeval


! ************************************************************************************ !
 
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!!					!!
!!		SPFitA0		!!	
!!					!!
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!

! Additional subroutine written as a preliminary step to SPFitA for Fortran 95 implementation.
! This ensures that all array dimensions are declared in SPFitA.

! Aubhik Khan and Julia K. Thomas 20.03.00

! Dimensions
!
!				m		number of dimensions of the multivariate function
!				rknots	vector of total knot points (interior knots + 2) in each dimension
!				clmat	clmat + 2 is the largest number of knots in any dimension
!				numelem	4 x m

SUBROUTINE SPFitA0(m, rknots, clmat, numelem)

INTEGER(ik):: m, rknots(m), rknotspoly(m), rknotsint(m), clmat, numelem(4,m), vars
INTENT(IN):: rknots, m
INTENT(OUT):: clmat, numelem

rknotspoly = rknots - 1_ik; rknotsint = rknots - 2_ik
clmat = MAXVAL(rknotsint)

! compute dimensions of matrices of spline coefficients
numelem = 4_ik
numelem(1,1) = PRODUCT(rknots)/rknots(1)

DO vars = 1_ik, m-1_ik, 1_ik
	numelem(2, vars) = rknots(vars)
	numelem(3,vars) = rknotspoly(vars)*numelem(1,vars)
	numelem(1, vars + 1) = (numelem(3, vars)*numelem(4, vars))/rknots(vars+1)
END DO

numelem(2, m) = rknots(m)
numelem(3,m) = rknotspoly(m)*numelem(1,m)

END SUBROUTINE SPFitA0

! ************************************************************************************ !
 
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!!					!!
!!		SPFitA		!!	
!!					!!
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!

! This subroutine generates a tridiagonal factorisation of the LHS matrix of coefficients 
! which are used to solve for the interior slopes of the univariate piecewise polynomial 
! splines.  As this procedure is required for each dimension of the multivariate spline, 
! it accomplishes this by undertaking repeated calls to SPLHS.  The relevant matrices for 
! each dimension are recorded in LMAT and UMAT.  Differences of the knots, used to generate
! divided differences, are also produced and, for each dimension, stored in DTAUMAT.  This 
! subroutine is required only once for any set of functions defined upon the same knot 
! product.

! Dimensions
!
!				LMAT	(m x clmat)
!				UMAT	(2 x (clmat+1) x m)
!				DTAUMAT (clmat+1 x numelem(1,m) x m)

SUBROUTINE SPFitA(m, rknots, clmat, knots, numelem, LMAT, UMAT, DTAUMAT)

INTEGER(ik):: m, rknotsint(m), rknotspoly(m), rknots(m), numelem(4,m), vars, &
	clmat, r, i
REAL(rk):: knots(clmat+2,m), LMAT(m, clmat), UMAT(2, clmat+1, m), DTAUMAT(clmat+1, numelem(1,m), m)
REAL(rk), DIMENSION(:,:), ALLOCATABLE:: L0, U0, dtau, dtau0
CHARACTER(20):: indicator

INTENT(IN):: rknots, m, knots, clmat, numelem
INTENT(OUT):: LMAT, UMAT, DTAUMAT

rknotspoly = rknots - 1_ik; rknotsint = rknots - 2_ik

! Not-a-knot endpoint condition
indicator = 'not-a'

LMAT = 0.0_rk
UMAT = 0.0_rk
DTAUMAT = 0.0_rk

! Allocate tridiagonal factorisation matrices for LHS coefficients for 
! spline slopes, L0, U0 and divided differences dtau for the first 
! dimension.
r = rknotsint(1)
ALLOCATE(L0(1, r-1), U0(2, r), dtau(r+1,1), dtau0(r+1,numelem(1,1)))

CALL SPLHS(knots(1:rknots(1),1), r, indicator, L0, U0, dtau)

FORALL (vars = 1:numelem(1,1)) dtau0(:,vars) = dtau(:,1)

LMAT(1,1:rknotsint(1) - 1) = L0(1,:)
UMAT(1:2, 1:rknotsint(1),1) = U0
DTAUMAT(1:rknotspoly(1), 1:numelem(1,1),1) = dtau0

! Repeat the above procedure for each additional dimension, recording all 
! matrices in LMAT, UMAT and DTAUMAT.
DEALLOCATE(L0, U0, dtau, dtau0)

DO i = 2, m, 1

	r = rknotsint(i)
	ALLOCATE(L0(1, r-1), U0(2, r), dtau(r+1,1), dtau0(r+1,numelem(1,i)))
	CALL SPLHS(knots(1:rknots(i),i), r, indicator, L0, U0, dtau)
	FORALL (vars = 1:numelem(1,i)) dtau0(:,vars) = dtau(:,1)
	LMAT(i,1:rknotsint(i) - 1) = L0(1,:)
	UMAT(1:2, 1:rknotsint(i),i) = U0(:,:)
	DTAUMAT(1:rknotspoly(i), 1:numelem(1,i),i) = dtau0
	DEALLOCATE(L0, U0, dtau, dtau0)

END DO

END SUBROUTINE SPFitA


! ************************************************************************************ !
 
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!!					!!
!!		SPFitB		!!	
!!					!!
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!

! This subroutine applies the piecewise cubic polynomial tensor product method to fitting
! multivariate splines to a function whose values are recorded in F0.  It accomplishes 
! this through repeated calls to SPpp.  

! Dimensions
!
!               m       scalar
!               c1mat   scalar
!               rknots  m vector
!               numelem 4 by m
!				LMAT	(m x clmat)
!				UMAT	(2 x (clmat+1) x m)
!				DTAUMAT (clmat+1 x numelem(1,m) x m)
!				F0		rknots(1) x [rknots(2)*...*rknots(m)] where rknots(i) -2 is the number of 
!															  number of interior knots in i-th 
!															  dimension, i = 1,...,m.  This is 
!                                                             simply rknots(1) x numelem(1,1).
!				ctemp	(4 x numelem(3,m))


SUBROUTINE SPFitB(m, c1mat, rknots, numelem, LMAT, UMAT, DTAUMAT, F0, ctemp)

INTEGER(ik):: m, c1mat, rknots(m), rknotspoly(m), rknotsint(m), r, i, gathernow, rearrange, numelem(4,m)
REAL(rk):: LMAT(m, c1mat), F0(rknots(1), numelem(1,1)), UMAT(2, c1mat+1, m), DTAUMAT(c1mat+1, numelem(1,m), m)
REAL(rk):: ctemp(4, numelem(3,m))
REAL(rk), ALLOCATABLE:: L0(:,:), Wkw(:,:), WK(:,:), ctemp0(:,:), SVEC(:,:), takeout(:,:), junk(:,:)
REAL(rk), ALLOCATABLE:: U0(:,:,:), dtau0(:,:,:) 
CHARACTER(20):: indicator

INTENT(IN):: LMAT, UMAT, DTAUMAT, numelem, rknots, m, F0
INTENT(OUT):: ctemp

rknotspoly = rknots - 1_ik; rknotsint = rknots - 2_ik

! Fit a univariate spline in the first dimension, corresponding to the rows of F
indicator = 'not-a'

r = rknotsint(1)

ALLOCATE(Wkw(rknots(1), numelem(1,1)), L0(1, r - 1), U0(2, r, 1), &
	dtau0(r+1, numelem(1,1),1), SVEC(2, numelem(1,1)), ctemp0(4, numelem(3,1)))

Wkw = F0
L0(1,:) = LMAT(1,1:r-1)
U0(:,:,1) = UMAT(1:2, 1:r,1)
dtau0(:,:,1) = DTAUMAT(1:r+1, 1:numelem(1,1),1)
SVEC = 0.0_rk
ctemp0 = 0.0_rk

CALL SPpp(Wkw, r, numelem(1,1), L0, U0, dtau0, indicator, SVEC, ctemp0) 

! Fit univariate splines in additional dimensions applying the tensor product approach
DO i = 2, m, 1

	gathernow = numelem(3, i-1)/numelem(2,i)

	ALLOCATE(takeout(4, gathernow), WK(rknots(i), 4*gathernow), junk(1,4*gathernow))

	DO rearrange = 1, rknots(i), 1
		
		takeout = ctemp0(:, 1 + (rearrange-1)*gathernow : rearrange*gathernow)
		junk = RESHAPE(takeout, (/1, 4*gathernow/))
		WK(rearrange, 1 : 4*gathernow) = junk(1,:)

	END DO

	DEALLOCATE(L0, U0, dtau0, SVEC, ctemp0, takeout, junk)

	r = rknotsint(i)

	ALLOCATE(L0(1, r - 1), U0(2, r, 1), dtau0(r+1, numelem(1,i),1), &
		SVEC(2, numelem(1,i)), ctemp0(4, numelem(3,i)))

	L0(1,:) = LMAT(i,1:r-1)
	U0(:,:,1) = UMAT(1:2, 1:r,i)
	dtau0(:,:,1) = DTAUMAT(1:r+1, 1:numelem(1,i),i)
	SVEC = 0.0_rk
	ctemp0 = 0.0_rk

	CALL SPpp(WK, r, numelem(1,i), L0, U0, dtau0, indicator, SVEC, ctemp0) 

	DEALLOCATE(WK)

END DO

ctemp = ctemp0

DEALLOCATE( L0, U0, dtau0, SVEC, ctemp0, Wkw)

END SUBROUTINE SPFitB

! ************************************************************************************ !
 
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!!					!!
!!	FastSplinEval	!!	
!!					!!
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!

! Dimensions
!
!	m			the number of dimensions
!	ccols		c2mat = numelem(3,m) - the number of columns in cwork, described as ctemp above
!	cwork		4 x ccols
!	rknotspoly	a vector of length m, rknots - 1
!   knotsrows	c1mat + 2, the largest number of knots in any dimension, the number of rows of knotsvec
!   knotsvec	knotsrows x m arrar, each column being knots in the i-th dimension, i = 1,...,m.
!   xvec		a vector of length m containing the point to be evaluated using the cwork spline defined
!				on the knotsvec multivariate knots
!	orders		a vector of length m with the order in each direction, orders(i) = 0,1, or 2, i = 1,..,m.

PURE FUNCTION FastSplinEval(m, ccols, cwork, rknotspoly, knotsrows, knotsvec, xvec, orders)

USE KindSet

IMPLICIT NONE

INTEGER:: coljtake, ccols, knotsrows, dennum
INTEGER(ik):: m, rknotspoly(m), i, j, coeffs(m), colnum, numberofcoeffs, &
	knothigh, knotlow, knotplace, s, knotloc(m), orders(m)
INTEGER(ik), DIMENSION(:,:), ALLOCATABLE:: jreach, jtake
REAL(rk):: cwork(4,ccols), knotsvec(knotsrows,m), xvec(m), varval, FastSplinEval
REAL(rk), DIMENSION(:,:), ALLOCATABLE:: intercoeffs
REAL(rk), DIMENSION(:), ALLOCATABLE:: ctransfer
INTENT(IN):: m, cwork, rknotspoly, knotsvec, xvec, orders, ccols, knotsrows

DO i = 1, m, 1

	knothigh = rknotspoly(i) + 1						! total number of knots = number of 											
	knotlow = 1											! polynomial pieces + 1
		
  DO 

  IF (knothigh - knotlow.LE.1_ik) THEN					! exit indefinite do loop (while loop)
	EXIT
  END IF      

     knotplace = (knothigh + knotlow)/2_ik				! bisect
	
     IF (knotsvec(knotplace,i).GE.xvec(i)) THEN			! is this tau(i) >= x(i)?
        knothigh = knotplace;							! yes, hence it is the upper bound on the bracket
     ELSE
        knotlow = knotplace;							! no, hence it is the lower bound on the bracket
     END IF
	      
  END DO													

	! the polynomial piece for x(i)
	knotloc(i) = knotlow

END DO

coeffs = 0_ik
dennum = MINVAL(rknotspoly)
colnum = size(cwork, 2)
coljtake = 4_ik**(m-1_ik)

ALLOCATE(jreach(m, colnum/dennum), jtake(m, coljtake))

jreach = 0_ik
jtake = 0_ik

DO i = m,1,-1
   
   coeffs(i) = colnum/rknotspoly(i)
   colnum = coeffs(i)/4
	
	DO j = 0,coeffs(i)-1, 1
		jreach(i,j+1) = j
	END DO

   jreach(i,1:coeffs(i)) = jreach(i,1:coeffs(i))*rknotspoly(i)
   jreach(i,1:coeffs(i)) = knotloc(i) + jreach(i,1:coeffs(i))

END DO

jtake(1,1) = jreach(1,coeffs(1))

DO i = 2,m,1
   
   numberofcoeffs = 4**(i-2)
   
   DO s = 1, numberofcoeffs, 1 
      
      jtake(i,(s-1)*4 + 1:s*4) = jreach(i,(jtake(i-1,s)-1)*4 + 1:(jtake(i-1,s)-1)*4 + 4)
      
   END DO
   
END DO

ALLOCATE(intercoeffs(4, coljtake))

intercoeffs = cwork(1:4, jtake(m,:))

DO i = m, 1, -1
   
   colnum = size(intercoeffs,2)

   ALLOCATE(ctransfer(colnum))

   ctransfer = 0.0_rk

   varval = xvec(i) - knotsvec(knotloc(i),i)
     
   IF (orders(i).EQ.0_ik) THEN

	ctransfer = intercoeffs(4,:)

      DO j = 3, 1, -1
         ctransfer = ctransfer*varval + intercoeffs(j,:)
      END DO
      
   ELSEIF (orders(i).EQ.1_ik) THEN
      
	  ctransfer = (3.0_rk*(varval**2.0_rk))*intercoeffs(4,:) &
		+ (2.0_rk*varval)*intercoeffs(3,:) &
		+ intercoeffs(2,:)
 
   ELSE
      
	  ctransfer = (6.0_rk*varval)*intercoeffs(4,:) + 2.0_rk*intercoeffs(3,:) 
     
   END IF    
   
	IF (colnum.GT.1_ik) THEN
	
		DEALLOCATE(intercoeffs)
		ALLOCATE(intercoeffs(4, colnum/4))
		colnum = colnum/4
		intercoeffs = RESHAPE(ctransfer, (/4, colnum/))
		DEALLOCATE(ctransfer)

	ELSE

		FastSplinEval = ctransfer(1)

	END IF
  
END DO

DEALLOCATE(ctransfer, intercoeffs, jreach, jtake)

END FUNCTION FastSplinEval

END MODULE PPSplineFit3edit