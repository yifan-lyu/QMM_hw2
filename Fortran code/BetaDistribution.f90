SUBROUTINE BetaDistribution(alpha0, beta0, xb, k, p0, expectation, probability)

! subroutine betafunction.for in FORTRAN77 from mincob taken from 
! http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html

USE KINDSET

IMPLICIT NONE

INTEGER(ik):: expectation
REAL(rk):: alpha0, beta0, xb, k, p0, integral, probability, newbeta, orgbeta

INTENT(IN):: alpha0, beta0, xb, k, p0, expectation
INTENT(OUT):: probability

EXTERNAL INCOBF, BETAF, GAMMAF

probability = 0.0_rk

IF (expectation.EQ.0_ik) THEN								! Return Pr(0.LT.x.LT.xb)


	IF (xb.EQ.0.0_rk) THEN
		
		probability = p0									! possible mass point at lower bound of 0

	ELSE

		CALL INCOBF(alpha0, beta0, xb/k, integral)			! incomplete beta function scaled for k
		probability = p0 + (1.0_rk - p0)*integral			! cumulativbe density given xb upper bound

	END IF


ELSE														! Return conditional expectation


	CALL BETAF(alpha0, beta0, newbeta)						! Adjust for denominator term in incomplete 
	CALL BETAF(alpha0+1.0_rk, beta0, orgbeta)				! Beta function.
	CALL INCOBF(alpha0+1.0_rk, beta0, xb/k, integral)	    ! Exploit alpha1 = alpha0+1 would give expectation  
															! if denominator complete BETA were not adjusted.			
	probability = integral*orgbeta/newbeta				
	probability = k*(1.0_rk - p0)*probability

END IF 

END SUBROUTINE BetaDistribution