PROGRAM BenchmarkPlanner

USE KINDSET
USE ppsplinefit3edit

! This is a simplified version of the program used to compute the benchmark economy for 
! Khan and Thomas (2007).  It solves the benchmark model without inventories using a 
! planners problem.  The recursive competitive equilibrium programs are not included.  
! 
! The program family is 
! BenchmarkPlanner.F90 to solve the model 
!	internal subroutines: Setup, Rouwenhorst
!   external subroutines: DR (solve), lambdabound, Linspace, Logspace, Tauchen2, 
!						  Normal, erf.f
! BenchmarkSim.F90 to simulate the model 
!   internal subroutines: Zshock, StdDev, Correlation
!   external subroutines: DR, HPFilterS.F90
! Modules: kindset.f90, ppslinefit3edit.f90

IMPLICIT NONE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                               !
!   Interface Blocks for External Subroutines   !   
!                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

interface

SUBROUTINE lambdabound(z0val, gz0val, kval, ktarget, alpha, Apref, delta, thetam, thetan, Atech, lambda, precision)

use kindset

REAL(rk):: z0val, gz0val, kval, ktarget, alpha, Apref, delta, thetam, thetan, Atech, &
           lambda, precision
                      
end subroutine lambdabound

           
SUBROUTINE DR(alpha, thetam, thetan, Atech, delta, Apref, rknotsk, &
			  knotsk, csk, zval, gzval, kval, lambda, musol, &
			  csol, msol, nmsol, nksol, kfsol, ysol, switch, erreqc, precision, lambdalow, &
			  lambdahigh)

use kindset

INTEGER(ik):: rknotsk 
REAL(rk):: alpha, thetam, thetan, Atech, delta, Apref, &
		   knotsk(rknotsk), csk(4,rknotsk-1), zval, gzval, kval, lambda, musol, csol, &
		   msol, nmsol, nksol, kfsol, precision, lambdalow, lambdahigh, &
		   erreqc, switch, ysol

end subroutine dr


SUBROUTINE BenchmarkSim(Nz, Z, PI, rknotsk, alpha, thetam, thetan, &
                        Atech, delta, beta, precisiond, Apref, knotsk, cskd, kbound, tauchenscale1, tauchenscale2, &
                        datafile, shockfile, simlength, again, ZIsim, Qsim, Psim, Ksim)
use kindset

INTEGER(ik):: Nz, rknotsk, simlength, ZIsim(simlength), again
REAL(rk):: Z(2, Nz), PI(Nz, Nz), cskd(4,rknotsk-1, Nz), &
		   alpha, thetam, thetan, Atech, delta, beta, precisiond, Apref, &
		   kbound(2), knotsk(rknotsk), Zsim(2,simlength), Ksim(simlength+1), Qsim(simlength), &
		   Psim(simlength), tauchenscale1, tauchenscale2
CHARACTER(30):: datafile, shockfile, resultfile		   

end subroutine benchmarksim		   

end interface		   



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                           !
!   Main Program BenchmarkPlanner.f90 Type Declarations     !
!                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		   
		   
INTEGER(ik):: Nz2(2), Nz, rknotsk, rknotskpoly, rknotskint, iz, ik0, iteration, simlength, again, time, &
              iz1, iz2, iz1f, iz2f
INTEGER(ik), ALLOCATABLE:: zIsim(:)
REAL(rk):: alpha, thetam, thetan, Atech, delta, rhoz2(2), stde2(2), tauchenscale1, &
		   tauchenscale2, beta, kbound(2), precisionv, precisiond, SVEC(2,1), Apref, zval, gzval, zmval, &
		   kval, kchoice(1), kfsol, csol, lambdasol, musol, nmsol, msol, nksol, epsilonv, errorval, q, &
		   lsol, switch, VX(3,1), klow, khigh, Vf, epsilond, ysol, begin, finish, p, DV0
REAL(rk), ALLOCATABLE:: Z(:,:), PI(:,:), knotsk(:), LW(:,:), UW(:,:), dtauW(:), csk(:,:, :), &
						V(:,:), TV(:,:), Kfvec(:,:), Lambdavec(:,:), EV(:), cskd(:,:,:), &
						DV(:,:), TDV(:,:), EDV(:), LL(:,:), LH(:,:), Wcaliph(:,:), QB(:), &
						PB(:), KB(:), csk0(:,:,:), Derror(:,:), D2VLook(:,:), &
						Z1vec(:), Z2vec(:), PI1vec(:,:), PI2vec(:,:)
CHARACTER(30):: datafile, resultfile, simulation
CHARACTER(20):: indicator
EXTERNAL LINSPACE, LOGSPACE, Tauchen2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							!
!	Setup of the Program	!
!							!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL CPU_TIME(begin) ! We time all calculations 

CALL SETUP(alpha, thetam, thetan, Atech, delta, rhoz2, stde2, Apref, beta, &
		   switch, precisionv, precisiond, Nz2, rknotsk, kbound, datafile, resultfile)

WRITE(*,*) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							        !
!	Discretize two shock processes  !
!							        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Two independent shocks to intermediate goods and final goods firms, Nz is total number of 
! exogenous states
Nz = Nz2(1)*Nz2(2)

ALLOCATE(Z1vec(Nz2(1)), PI1vec(Nz2(1), Nz2(1)), Z2vec(Nz2(2)), PI2vec(Nz2(2), Nz2(2)), &
         Z(2,Nz), PI(Nz,Nz))

! The benchmark economy in the paper sets tauchenscale = 2.0 for z1
WRITE(*,'(1X,A)', ADVANCE = 'NO') ' enter std z1 (z) scale factor for the Tauchen grid (typically 2.0 ): '
READ(*,*)  tauchenscale1

IF (tauchenscale1.NE.0.0_rk.AND.Nz2(1).GT.1_ik) THEN
	CALL TAUCHEN2(0.0_rk, stde2(1), rhoz2(1), tauchenscale1, Nz2(1), Z1vec, PI1vec)	
ELSE 
	CALL Rouwenhorst(rhoz2(1), stde2(1), Nz2(1), Z1vec, PI1vec)	
END IF

Z1vec = DEXP(Z1vec)

! The benchmark economy sets tauchenscale2 = 0.0 (no second shock)
WRITE(*,'(1X,A)', ADVANCE = 'NO') ' enter std z2 (gz) scale factor (typically 0.0 to use Rouwenhorst): '
READ(*,*)  tauchenscale2

IF (tauchenscale2.NE.0.0_rk.AND.Nz2(2).GT.1_ik) THEN
    CALL TAUCHEN2(0.0_rk, stde2(2), rhoz2(2), tauchenscale2, Nz2(2), Z2vec, PI2vec)
ELSE
    CALL Rouwenhorst(rhoz2(2), stde2(2), Nz2(2), Z2vec, PI2vec)
END IF

Z2vec = DEXP(Z2vec)

! Place the first shock, z, in the first row of Z, the second shock gz in the second row.
DO iz1 = 1, Nz2(1)
    DO iz2 = 1, Nz2(2)
        Z(1,(iz1-1)*Nz2(2) + iz2) = Z1vec(iz1)
        Z(2,(iz1-1)*Nz2(2) + iz2) = Z2vec(iz2)
    END DO
END DO

DO iz1 = 1, Nz2(1)
    DO iz2 = 1, Nz2(2)
        DO iz1f = 1, Nz2(1)
            DO iz2f = 1, Nz2(2)
                PI((iz1 -1)*Nz2(2) + iz2, (iz1f-1)*Nz2(2) + iz2f)  = PI1vec(iz1,iz1f)*PI2vec(iz2, iz2f)
            END DO      
        END DO
    END DO
END DO        
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !    
!   Univariate piecewise polynomial cubic spline interpolation      !
!   for the Value function and its derivative with respect to       !
!   capital.  See notes on spline interpolation                     !  
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! UNIVARIATE KNOT SETUP.  The dimension of LW, UW, dtauW and csk should be robust 
! for any one-dimensional example.
rknotskpoly = rknotsk - 1_ik; rknotskint = rknotsk - 2_ik

ALLOCATE(knotsk(rknotsk), LW(1, rknotskint-1), UW(2, rknotskint), dtauW(rknotskpoly))

CALL LINSPACE(kbound(1), kbound(2), rknotsk, knotsk)

! Endpoint conditions for univariate spline approximation.
indicator = '                   '
indicator(1:10) = 'not-a-knot'	! endpoint selection for univariate spline coefficients.
SVEC = 0.0_rk

CALL SPLHS(knotsk, rknotskint, indicator, LW, UW, dtauW)

ALLOCATE(csk(4, rknotskpoly, Nz), cskd(4, rknotskpoly, Nz), csk0(4, rknotskpoly, Nz))

! spline coefficents for the value function and its derivative will be held in csk and cskd.
csk = 0.0_rk; cskd = 0.0_rk

! Array allocations for contraction mapping
ALLOCATE(V(Nz, rknotsk), TV(Nz, rknotsk), Kfvec(Nz, rknotsk), Lambdavec(Nz, rknotsk), &
		 EV(rknotsk), DV(Nz, rknotsk), TDV(Nz, rknotsk), EDV(rknotsk), LL(Nz, rknotsk), &
		 LH(Nz, rknotsk), Wcaliph(Nz, rknotsk), Derror(Nz, rknotsk), D2VLook(Nz, rknotsk))
		 
! Initial value function
DO iz = 1, Nz
	zval = Z(1,iz)  ! Take the first shock for z here
	DO ik0 = 1, rknotsk
		kval = knotsk(ik0)
		V(iz, ik0) = DLOG(zval*(kval**alpha))/(1.0_rk - beta)
		DV(iz,ik0) = (alpha/(1.0_rk - beta))*(kval**(alpha - 1.0_rk))
	END DO
END DO

epsilonv = 2.0_rk*precisionv; iteration = 0_ik


!!!!!!!!!!!!!!!!!!!!!
!					!
!	lambda bounds	! 
!					!
!!!!!!!!!!!!!!!!!!!!!

! Determine low and high values of lambda, the marginal utility of current consumption, so 
! that for each element of the aggregate state (z(1), z(2), k), the choice of k(t+1) does
! not violate the knots on k.  The values of lambda are implied by the first-order condition
! for k(t+1).    
klow = knotsk(1)
khigh = knotsk(rknotsk)
DO iz = 1,Nz
	zval = Z(1,iz); gzval = Z(2,iz)
	DO ik0 = 1, rknotsk
		kval = knotsk(ik0)	
		CALL lambdabound(zval, gzval, kval, klow, alpha, Apref, delta, thetam, thetan, Atech, LL(iz,ik0), precisiond)
		CALL lambdabound(zval, gzval, kval, khigh, alpha, Apref, delta, thetam, thetan, Atech, LH(iz,ik0), precisiond)	
	END DO
END DO



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							!
!	Contraction mapping		!
!							!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO

	iteration = iteration + 1_ik

	IF (epsilonv.LT.precisionv) THEN
		EXIT
	END IF
	
	! Fit spline to expected value function, one for each current z 
	! that yields W(k) = E[V(z',k')|z], where W(k) = V(k,z).
	DO iz = 1, Nz		
		DO ik0 = 1, rknotsk
			EV(ik0) = DOT_PRODUCT(PI(iz,:), beta*V(:,ik0))
			EDV(ik0) = DOT_PRODUCT(PI(iz,:), beta*DV(:,ik0))
		END DO
	
		CALL SPpp(EV, rknotskint, 1, LW, UW, dtauW, indicator, (/-1.0_rk, -0.5_rk/), csk(:,:,iz))
		CALL SPpp(EDV, rknotskint, 1, LW, UW, dtauW, indicator, (/0.0_rk, 0.0_rk/), cskd(:,:,iz))		
	END DO
		
    ! Solve the optimisation problem
	DO iz = 1, Nz
		zval = Z(1,iz); gzval = Z(2,iz); zmval = zval
		DO ik0 = 1, rknotsk
			kval = knotsk(ik0); errorval = 0.0_rk

			CALL DR(alpha, thetam, thetan, Atech, delta, Apref, rknotsk, &
					knotsk, cskd(:,:,iz), zval, gzval, kval, lambdasol, musol, &
					csol, msol, nmsol, nksol, kfsol, ysol, switch, errorval, precisiond, &
					LL(iz,ik0), LH(iz,ik0))
		
			lsol = 1.0_rk - (nmsol + nksol)
			
			! Second-order Taylor expansion to evaluate E[V(kf)|(z(1),z(2)] if kf 
			! violates the knots for k.  This is entirely unnecessary if lambdabounds
			! were used properly above.
			IF (kfsol.LT.klow) THEN	! 2nd order Taylor expansions if knotsk violations below
				kchoice(1) = klow
				CALL SPeval(csk(:,:,iz), knotsk, rknotskint, 1, kchoice, 1, 2, VX)
				Vf = VX(1,1) + VX(2,1)*(kfsol-klow) + switch*0.5_rk*VX(3,1)*(kfsol-klow)**2.0_rk
			ELSEIF(kfsol.GT.khigh) THEN ! or if knotsk violated above. 
				kchoice(1) = khigh	
				CALL SPeval(csk(:,:,iz), knotsk, rknotskint, 1, kchoice, 1, 2, VX)
				Vf = VX(1,1) + VX(2,1)*(kfsol-khigh) + switch*0.5_rk*VX(3,1)*(kfsol-khigh)**2.0_rk
			ELSE				 ! No knots violated.	
				kchoice(1) = kfsol
				CALL SPeval(csk(:,:,iz), knotsk, rknotskint, 1, kchoice, 1, 0, VX)
				Vf = VX(1,1)
			END IF			
			
			IF (errorval.GE.0.01_rk) THEN

			WRITE(*, FMT = '(1X, I3, A, F6.3, A, F6.3, A, F6.3, A,F6.3,A,F6.3, A, 3(F8.4), A, F8.4)', ADVANCE = 'YES') &
			iteration, ' (', zval, ',', gzval, ',', kval, ') ll, lh = ',LL(iz,ik0), ',', LH(iz,ik0),' (lambda, c, kf) = ', &
			lambdasol, csol, kfsol, ' error: ', errorval

			END IF

			TV(iz,ik0) = DLOG(csol) + Apref*(lsol) + Vf
			TDV(iz,ik0) = lambdasol*zmval*alpha*(kval**(alpha - 1.0_rk))*(nksol**(1.0_rk - alpha)) + musol*(1.0_rk - delta)
		
			p = 1.0_rk/csol; q = lambdasol/musol
		
			Kfvec(iz,ik0) = kfsol; Lambdavec(iz,ik0) = lambdasol

		END DO

	END DO

	epsilonv = MAXVAL(DABS(TV - V))
	epsilond = MAXVAL(DABS(TDV - DV))
	WRITE(*,'(1X,A,I4,A,2(E18.4))', ADVANCE = 'YES') ' SVEC used iteration: ', iteration, '     norm (V,DV) = ', epsilonv, epsilond
	V = TV; DV = TDV
		
END DO

CALL CPU_TIME(finish); finish = finish - begin


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								!
!	 Consistency of V and DV	!
!								!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Verify consistency between DV
WRITE(*,*) ' ' 
WRITE(*,*) ' Consistency between V and monotone operator on DV, followed by D2V(k,z) (each pair is for a z value) ' 

csk0 = 0.0_rk

DO iz = 1, Nz
	
	DO ik0 = 1, rknotsk
		EV(ik0) = V(iz,ik0)		
	END DO

	CALL SPpp(EV, rknotskint, 1, LW, UW, dtauW, indicator, SVEC, csk0(:,:,iz))

	DO ik0 = 1, rknotsk

		kval = knotsk(ik0)

		! Adapted for kchoice for CVF 6.6B 
		IF (kval.LT.klow) THEN	! 2nd order Taylor expansions if knotsk violations below
			kchoice(1) = klow
			CALL SPeval(csk0(:,:,iz), knotsk, rknotskint, 1, kchoice, 1, 2, VX)
			DV0 = VX(2,1) + VX(3,1)*(kval-klow)
		ELSEIF(kval.GT.khigh) THEN ! or if knotsk violated above. 	
			kchoice(1) = khigh
			CALL SPeval(csk0(:,:,iz), knotsk, rknotskint, 1, kchoice, 1, 2, VX)
			DV0 = VX(2,1) + VX(3,1)*(kval-khigh)
		ELSE				 ! No knots violated.	
			kchoice(1) = kval
			CALL SPeval(csk0(:,:,iz), knotsk, rknotskint, 1, kchoice, 1, 2, VX)
			DV0 = VX(2,1)
		END IF		

		Derror(iz,ik0) = DV0 - DV(iz,ik0)
		D2VLook(iz,ik0) = VX(3,1)

	END DO

END DO

DO ik0 = 1, rknotsk
	DO iz = 1, Nz
		WRITE(*,'(1X, F8.4)', ADVANCE = 'NO') Derror(iz,ik0) 
	END DO

	DO iz = 1, Nz
		WRITE(*,'(1X, F12.4)', ADVANCE = 'NO') D2VLook(iz,ik0) 
	END DO

		WRITE(*,*) ' '
END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							!
!	Planner Simulation		!
!							!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,'(1X,A,F6.3,A)', ADVANCE = 'NO') ' Value Function Convergence in ', finish, ' seconds. Simulate the economy? '
READ(*,*) simulation

IF (simulation.EQ.'y') THEN

	WRITE(*, '(1X,A)', ADVANCE = 'NO') ' simulation length (-1 to load existing simulation): '
	READ(*,*) simlength

	IF (simlength.EQ.-1) THEN
		simulation = 'n'	! Not a new simulation
		again = simlength
		WRITE(*, '(1X,A)', ADVANCE = 'NO') ' I have to know the actual existing simulation length: '
		READ(*,*) simlength
	END IF
	
	ALLOCATE(zIsim(simlength), QB(simlength), PB(simlength), KB(simlength+1))

	CALL BenchmarkSim(Nz, Z, PI, rknotsk, alpha, thetam, thetan, &
					 Atech, delta, beta, precisiond, Apref, knotsk, cskd, kbound, tauchenscale1, &
					 tauchenscale2, datafile, resultfile, simlength, again, zIsim, QB, PB, KB)

	
	IF (simulation.EQ.'y') THEN		! record a new simulation

		OPEN(UNIT = 30, FILE = resultfile, ACTION = 'READWRITE', STATUS = 'REPLACE')	

		! Record the results as an initial condition for the dynamic inventory model
	
		WRITE(30,*) simlength
	
		DO time = 1,simlength
			WRITE(30,*) zIsim(time), QB(time), PB(time)
		END DO

		DO time = 1, simlength + 1
			WRITE(30,*) KB(time)
		END DO

		WRITE(30,*) tauchenscale1   ! spurious write statements for backward consistency
		WRITE(30,*) tauchenscale1
		WRITE(30,*) tauchenscale1
		WRITE(30,*) Nz
		
		DO iz1 = 1, Nz
		    WRITE(30,*) Z(1,iz1), Z(2,iz1)
		END DO
		
		DO iz1 = 1, Nz
		    DO iz2 = 1, Nz
		        WRITE(30,*) PI(iz1, iz2)
		    END DO
		END DO
		
		CLOSE(30)

	ELSE

		WRITE(*,*) ' Did not overwrite existing simulation file '

	END IF


END IF

 WRITE(*,'(1X,A)', ADVANCE = 'NO') ' The Baseline Planner New programs using ppslinefit3edit have completed.  Enter to exit. '
 READ(*,*)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!							!!
!!							!!
!!	 Internal Subroutines	!!
!!							!!
!!							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

!!!!!!!!!!!!!
!			!
!	SETUP	!
!			!
!!!!!!!!!!!!!

SUBROUTINE SETUP(alpha, thetam, thetan, Atech, delta, rhoz, stde, Apref, &
				beta, switch, precisionv, precisiond, Nz, Nk, kbound, datafile, resultfile)  

INTEGER(ik):: Nz(2), Nk, i, what
REAL(rk):: alpha, thetam, thetan, Atech, delta, rhoz(2), stde(2), &
		   beta, kbound(2), precisionv, precisiond, Apref, switch
CHARACTER(30):: datafile, resultfile

INTENT(OUT):: alpha, thetam, thetan, Atech, delta, rhoz, stde, beta, switch, &
			  precisionv, precisiond, Nz, Nk, kbound, datafile, resultfile, Apref

DO 

	WRITE(*,FMT = '(1X,A)', ADVANCE = 'NO') ' enter parameter file: '
	READ(*,*) datafile
	
	
	OPEN (UNIT = 57, FILE = datafile, STATUS = 'OLD', ACTION = 'READ', IOSTAT = what, ERR = 100)
	
	! The subroutine goes to line 100, right below, if there's a file opening error.  
	100	IF (what.NE.0_ik) THEN 
			WRITE(*,*) ' Hey, I think I could not locate ', datafile
			WRITE(*,*) ' I got a file I/O error ', what
		ELSE
			EXIT
		END IF

END DO


OPEN (UNIT = 57, FILE = datafile, STATUS = 'OLD', ACTION = 'READ')

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
READ(57,*) (rhoz(i), i = 1,2)               			! peristence of z, z' = rhoz*z + epsilon	
READ(57,*)
READ(57,*) (stde(i), i = 1,2)							! variance of epsilon above	
READ(57,*)
READ(57,*) beta											! discount factor
READ(57,*)
READ(57,*) Apref										! preference term for leisure
READ(57,*)
READ(57,*) (kbound(i), i = 1,2)							! knot endpoints on capital
READ(57,*)
READ(57,*) (Nz(i), i = 1,2)								! number of total z values for each shock
READ(57,*)
READ(57,*)  Nk											! total number of knots, including endpoints, on k	
READ(57,*)
READ(57,*)  switch										! set to 1 for second order Taylor expansions
READ(57,*)  
READ(57,*)  precisionv									! precision term for value function convergence
READ(57,*)  
READ(57,*)  precisiond									! precision term for consumption determination
READ(57,*)  
READ(57,*)  resultfile

CLOSE (57)						! close datafile

END SUBROUTINE SETUP

!!!!!!!!!!!!!!!!!!!!!
!					!	
!					!
!	ROUWENHORST		!
!					!
!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Rouwenhorst(rho, sigmas, N, Z, PI)

! This subroutine is taken from AR1MARKOV.F90 and Rouwenhorst.m, 

INTEGER(ik), INTENT(IN):: N
INTEGER(ik)::i
REAL(rk):: p, q, zvar, epsilon, Y(N)
REAL(rk), INTENT(IN):: rho, sigmas
REAL(rk), INTENT(OUT):: Z(N), PI(N,N)
REAL(rk), DIMENSION(:,:), ALLOCATABLE:: hlag, h

! Retireive p and q using algorithm in Rouwenhorst.m
p = (rho + 1.0_rk)/2.0_rk
q = p

ALLOCATE(hlag(1,1))

hlag(1,1) = 1.0_rk

DO i = 2,N

	ALLOCATE (h(i,i))

	h = 0.0_rk
	h(1:i-1,1:i-1) = p*hlag
	h(1:i-1,2:i) = h(1:i-1,2:i) + (1.0_rk-p)*hlag
	h(2:i,1:i-1) = h(2:i,1:i-1) + (1.0_rk-q)*hlag
	h(2:i,2:i) = h(2:i,2:i) + q*hlag

	! Note that the following is valid despite i-1<2 
	! when i = 2 and h not being conformable with the 
	! scalar 2 at any time.
	h(2:i-1,:) = h(2:i-1,:)/2.0_rk

	DEALLOCATE(hlag)
	
	ALLOCATE(hlag(i,i))
	
	hlag = h
	
	DEALLOCATE(h)

END DO

PI = hlag

DEALLOCATE(hlag)	

zvar = (sigmas**2.0_rk)/(1.0_rk - rho**2.0_rk)
epsilon = DSQRT((N-1.0_rk)*zvar)

! These lines are adapted from Linspace.f90
DO i = 1_ik, N - 2_ik
	Y(i) = dble(i )
END DO

Y = ((2.0_rk*epsilon)/dble(N-1_ik))*Y - epsilon 

Z(1_ik) = -1.0_rk*epsilon
Z(2_ik:N-1_ik) = Y
Z(N) = epsilon

END SUBROUTINE Rouwenhorst 

END PROGRAM BenchmarkPlanner



! ***************************************************** !
! ***************************************************** !
!														!
!														!
!														!		
!		External Subroutine for Decision Rules			!
!														!
!														!
!														!
! ***************************************************** !
! ***************************************************** !


SUBROUTINE DR(alpha, thetam, thetan, Atech, delta, Apref, rknotsk, &
			  knotsk, csk, zval, gzval, kval, lambda, musol, &
			  csol, msol, nmsol, nksol, kfsol, ysol, switch, erreqc, precision, lambdalow, &
			  lambdahigh)

USE KINDSET
USE ppsplinefit3edit

IMPLICIT NONE

INTEGER(ik):: rknotsk, rknotskpoly, rknotskint, s1
REAL(rk):: alpha, thetam, thetan, Atech, delta, Apref, &
		   knotsk(rknotsk), csk(4,rknotsk-1), zval, gzval, kval, lambda, musol, csol, &
		   msol, nmsol, nksol, kfsol, precision, lambdalow, lambdahigh, r, kchoice(1), &
		   a0, b0, a, b, c, d, gzgsc, klow, khigh, hc(3,1), hd(3,1), erreqc, switch, fc, fd, &
		   eqclambda, ysol

INTENT(IN):: rknotsk, alpha, thetam, thetan, Atech, &
             delta, knotsk, csk, zval, gzval, kval, precision, switch, lambdalow, lambdahigh, Apref

INTENT(OUT):: musol, lambda, csol, msol, nmsol, nksol, kfsol, ysol, erreqc

rknotskint = rknotsk - 2_ik; rknotskpoly = rknotsk - 1_ik

! Upper and lower bounds for lambda, note mu now serves the role of the previous lambda.
klow = knotsk(1); khigh = knotsk(rknotsk)

! Golden section search on F.O.C. for k'
a0 = lambdalow; b0 = lambdahigh

a = a0; b = b0; erreqc = 97.0_rk
r = (3.0_rk - DSQRT(5.0_rk))/2.0_rk

! The c evaluation: F(c), this is now a little different, we get one mu from solve, and compare it 
! to the derivative of the value function with respect to capital.
c = a + r*(b-a); 

CALL solve(alpha, thetam, c, Atech, thetan, delta, Apref, zval, gzval, kval, msol, nmsol, &
		   nksol, musol, kfsol, csol, ysol)

IF (kfsol.LT.klow) THEN	! 2nd order Taylor expansions if knotsk violations below
	kchoice(1) = klow
	CALL SPeval(csk, knotsk, rknotskint, 1, kchoice, 1, 2, hc) 
	eqclambda = hc(1,1) + hc(2,1)*(kfsol - klow) + switch*0.5_rk*hc(3,1)*((kfsol - klow)**2.0_rk)
ELSEIF(Kfsol.GT.khigh) THEN ! or if knotsk violated above. 	
	kchoice(1) = khigh
	CALL SPeval(csk, knotsk, rknotskint, 1, kchoice, 1, 2, hc)
	eqclambda = hc(1,1) + hc(2,1)*(kfsol - khigh) + switch*0.5_rk*hc(3,1)*((kfsol - khigh)**2.0_rk)
ELSE				 ! No knots violated.	
	kchoice(1) = kfsol
	CALL SPeval(csk, knotsk, rknotskint, 1, kchoice, 1, 2, hc); eqclambda = hc(1,1)
END IF
fc = (eqclambda - musol)**2.0_rk

! The d evaluation
d = a + (1.0_rk - r)*(b-a);
CALL solve(alpha, thetam, d, Atech, thetan, delta, Apref, zval, gzval, kval, msol, nmsol, &
		   nksol, musol, kfsol, csol, ysol)
IF (kfsol.LT.klow) THEN	! 2nd order Taylor expansions if knotsk violations below
	kchoice(1) = klow
	CALL SPeval(csk, knotsk, rknotskint, 1, kchoice, 1, 2, hd) 
	eqclambda = hd(1,1) + hd(2,1)*(kfsol - klow)+ switch*0.5_rk*hd(3,1)*((kfsol - klow)**2.0_rk)
ELSEIF(kfsol.GT.khigh) THEN ! or if knotsk violated above. 	
	kchoice(1) = khigh
	CALL SPeval(csk, knotsk, rknotskint, 1, kchoice, 1, 2, hd) 
	eqclambda = hd(1,1) + hd(2,1)*(kfsol - khigh)+ switch*0.5_rk*hd(3,1)*((kfsol - khigh)**2.0_rk)
ELSE				 ! No knots violated.	
	kchoice(1) = kfsol
	CALL SPeval(csk, knotsk, rknotskint, 1, kchoice, 1, 1, hd); eqclambda = hd(1,1)
END IF
fd = (eqclambda - musol)**2.0_rk

DO
		  
	IF (DABS(d-c).LT.precision) THEN
		EXIT
	END IF
			                       
    s1 = s1 + 1_ik
               
    IF (fc.GE.fd) THEN
           
		gzgsc = c + (1.0_rk-r)*(b-c)
        ! [a c d b] <--- [c d z b]
        a = c
        c = d; fc = fd
		d = gzgsc 
		
		CALL solve(alpha, thetam, d, Atech, thetan, delta, Apref, zval, gzval, kval, msol, nmsol, &
		           nksol, musol, kfsol, csol, ysol)

		IF (kfsol.LT.klow) THEN	! 2nd order Taylor expansions if knotsk violations below
			kchoice(1) = klow
			CALL SPeval(csk, knotsk, rknotskint, 1, kchoice, 1, 2, hd) 
			eqclambda = hd(1,1) + hd(2,1)*(kfsol - klow) + switch*0.5_rk*hd(3,1)*((kfsol - klow)**2.0_rk)
		ELSEIF(Kfsol.GT.khigh) THEN ! or if knotsk violated above. 	
			kchoice(1) = khigh
			CALL SPeval(csk, knotsk, rknotskint, 1, kchoice, 1, 2, hd) 
			eqclambda = hd(1,1) + hd(2,1)*(kfsol - khigh) + switch*0.5_rk*hd(3,1)*((kfsol - khigh)**2.0_rk)
		ELSE				 ! No knots violated.	
			kchoice(1) = kfsol
			CALL SPeval(csk, knotsk, rknotskint, 1, kchoice, 1, 1, hd); eqclambda = hd(1,1)
		END IF
		fd = (eqclambda - musol)**2.0_rk; erreqc = fd; lambda = d
                    
	ELSE
		
		gzgsc = a + r*(d-a)
        ! [a c d b] <--- [a z c d]
        b = d
        d = c; fd = fc
        c = gzgsc 
		
		CALL solve(alpha, thetam, c, Atech, thetan, delta, Apref, zval, gzval, kval, msol, nmsol, &
		           nksol, musol, kfsol, csol, ysol)
		IF (kfsol.LT.klow) THEN	! 2nd order Taylor expansions if knotsk violations below
			kchoice(1) = klow
			CALL SPeval(csk, knotsk, rknotskint, 1, kchoice, 1, 2, hc) 
			eqclambda = hc(1,1) + hc(2,1)*(kfsol - klow) + switch*0.5_rk*hc(3,1)*((kfsol - klow)**2.0_rk)
		ELSEIF(Kfsol.GT.khigh) THEN ! or if knotsk violated above. 	
			kchoice(1) = khigh
			CALL SPeval(csk, knotsk, rknotskint, 1, kchoice, 1, 2, hc) 
			eqclambda = hc(1,1) + hc(2,1)*(kfsol - khigh) + switch*0.5_rk*hc(3,1)*((kfsol - khigh)**2.0_rk)
		ELSE				 ! No knots violated.	
			kchoice(1) = kfsol
			CALL SPeval(csk, knotsk, rknotskint, 1, kchoice, 1, 1, hc); eqclambda = hc(1,1)
		END IF
		fc = (eqclambda - musol)**2.0_rk; erreqc = fc; lambda = c            

    END IF
     
 
END DO


CONTAINS

! Determine decision rules given some value of lambda

SUBROUTINE solve(alpha, thetam, lambda, Atech, thetan,  delta, Apref, z0val, gz0val, kval, m, nm, nk, &
				 musol, kf, csol, y)
IMPLICIT NONE

REAL(rk):: alpha, thetam, lambda, Atech, thetan, Apref, z0val, gz0val, kval, m, nm, nk, &
		   musol, kf, csol, i, delta, g, y, D1R, zval 
INTENT(IN):: alpha, thetam, lambda, Atech, thetan, delta, Apref, z0val, gz0val, kval
INTENT(OUT):: m, nm, nk, musol, kf, csol, y

zval = z0val
g = gz0val; g = Atech*g

! Steps to obtain k'(lambda): (z,k) given, fix lambda, determine nk, m, nm
nk = (lambda*(1.0_rk - alpha)*zval/Apref)**(1.0_rk/alpha)
nk = nk*kval

m = zval*(kval**alpha)*(nk**(1.0_rk - alpha))

nm = zval*(1.0_rk - alpha)*(kval**alpha)*(nk**(-1.0_rk*alpha))
nm = (thetan/thetam)/nm
nm = nm*m

D1R = thetam*(m**(thetam - 1.0_rk))*(nm**thetan)	
y = g*(m**thetam)*(nm**thetan)

musol = lambda/(g*D1R)
csol = 1.0_rk/musol
i = y - csol
kf = (1.0_rk - delta)*kval + i

END SUBROUTINE solve

END SUBROUTINE DR


!!!!!!!!!!!!!!!!!!!!!
!					!
!	lambdabound		!
!					!
!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE lambdabound(z0val, gz0val, kval, ktarget, alpha, Apref, delta, thetam, thetan, Atech, lambda, precision)

USE KINDSET

IMPLICIT NONE

REAL(rk):: z0val, gz0val, kval, ktarget, alpha, Apref, delta, thetam, thetan, Atech, &
           lambda, kf, llow, lhigh, distance, flow, fhigh, f, precision, gz, zval
INTENT(IN):: z0val, gz0val, kval, ktarget, alpha, Apref, delta, thetam, thetan, Atech, &
             precision
INTENT(OUT):: lambda

zval = z0val
gz = gz0val; gz = Atech*gz

llow = precision
lambda = llow
CALL kfgivenlambda(alpha, gz, thetan, thetam, delta, Apref, zval, kval, lambda, kf)

flow = kf - ktarget

lhigh = 1.0_rk/precision
lambda = lhigh
CALL kfgivenlambda(alpha, gz, thetan, thetam, delta, Apref, zval, kval, lambda, kf)

fhigh = kf - ktarget

IF (flow*fhigh.GE.0.0_rk) THEN
	WRITE(*,'(1X,A,3(F8.4), A, 2(E12.4))', ADVANCE = 'YES')  ' CANNOT BISECT AT zval, kval, ktarget = ', zval, kval, &
	ktarget, ' flow, fhigh = ', flow, fhigh
	PAUSE
END IF

DO
	distance = lhigh - llow
	IF (distance.LT.precision) THEN
		EXIT
	END IF

	lambda = (llow + lhigh)/2.0_rk
	CALL kfgivenlambda(alpha, gz, thetan, thetam, delta, Apref, zval, kval, lambda, kf)
	f = kf - ktarget

	IF (f.GT.0.0_rk) THEN
		lhigh = lambda
	ELSE
		llow = lambda
	END IF

END DO


CONTAINS

SUBROUTINE kfgivenlambda(alpha, gz, thetan, thetam, delta, Apref, zval, kval, lambda, kf)

IMPLICIT NONE

REAL(rk):: alpha, gz, thetan, thetam, delta, Apref, zval, kval, lambda, kf, nk, m, &
		   nm, D1R, mu
INTENT(IN):: alpha, gz, thetan, thetam, delta, Apref, zval, kval, lambda
INTENT(OUT):: kf

! Steps to obtain k'(lambda): (z,k) given, fix lambda, determine nk, m, nm
nk = (lambda*(1.0_rk - alpha)*zval/Apref)**(1.0_rk/alpha); nk = nk*kval
m = zval*(kval**alpha)*(nk**(1.0_rk - alpha))
nm = zval*(1.0_rk - alpha)*(kval**alpha)*(nk**(-1.0_rk*alpha))
nm = (thetan/thetam)/nm
nm = nm*m

D1R = thetam*(m**(thetam - 1.0_rk))*(nm**thetan)	
mu = lambda/(gz*D1R)
kf = gz*(m**thetam)*(nm**thetan) + (1.0_rk - delta)*kval - 1.0_rk/mu

! End of Steps to obtain k'(lambda)

END SUBROUTINE kfgivenlambda

END SUBROUTINE lambdabound


!!!!!!!!!!!!!!!!!
!				!
!	LOGSPACE	!
!				!
!!!!!!!!!!!!!!!!!

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

! This is lifted directly from the linspace subroutine. 
DO j1 = 1, n - 2
	Y(j1) = DBLE(j1)
END DO

Y = ((ub - lb)/dble(n-1))*Y + lb

X(1) = lb
X(2:n-1) = Y
X(n) = ub

X = 10.0_rk**X

END SUBROUTINE LogSpace


!!!!!!!!!!!!!!!!!
!				!			
!	LinSpace	!
!				!
!!!!!!!!!!!!!!!!!

SUBROUTINE LinSpace(lb, ub, gridnum, X)

USE KINDSET

IMPLICIT NONE

INTEGER(ik):: gridnum, j1
REAL(rk):: lb, ub, X(gridnum), Y(gridnum-2_ik)

INTENT(IN):: lb, ub, gridnum
INTENT(OUT):: X

DO j1 = 1_ik, gridnum - 2_ik
	Y(j1) = dble(j1 )
END DO

Y = ((ub - lb)/dble(gridnum-1_ik))*Y + lb

X(1_ik) = lb
X(2_ik:gridnum-1_ik) = Y
X(gridnum) = ub

END SUBROUTINE LinSpace


!!	**************************************************************	!!
!!	**************************************************************	!!
!!																	!!	
!!																	!!
!!																	!!
!!																	!!
!!					DISCRETE APPROXIMATION OF A						!!
!!				 CONTINUOUS AUTOREGRESSIVE PROCESS					!!
!!																	!!
!!																	!!
!!																	!!
!!																	!!
!!	**************************************************************	!!
!!	**************************************************************	!!

! G. Tauchen (1986) 'Finite State Markov-Chain Approximations to 
! Univariate and Vector Autoregressions' Economics Leters 20: 177-181 

!  version 2 uses external normal function with erf.f fortran 77 subroutine 

SUBROUTINE TAUCHEN2(meaninnov, stdinnov, persistence, multiple, znum, ZVEC, PI)

USE KINDSET

IMPLICIT NONE

INTEGER(ik):: znum, j, gridsize, k
REAL(rk)::meaninnov, stdinnov, persistence, multiple, ZVEC(znum), PI(znum, znum), &
		  F1, F0, stdz, zlow, zhigh, meanz, w, z, lowerbound, Normal
EXTERNAL LINSPACE, Normal

INTENT(IN):: meaninnov, stdinnov, persistence, multiple, znum
INTENT(OUT):: ZVEC, PI

! set endpoints of stochastic variable grid at multiple m of the 
! standard deviation.  

stdz = stdinnov**2.0_rk
stdz = stdz/(1.0_rk - persistence**2.0_rk)
stdz = DSQRT(stdz)
meanz = meaninnov/(1.0_rk - persistence)
zlow = meanz - stdz*multiple
zhigh = meanz + stdz*multiple

lowerbound = meaninnov - stdinnov*DMAX1(10.0_rk, 2.0_rk*multiple)
gridsize = 10000

CALL LINSPACE(zlow, zhigh, znum, ZVEC)

PI = 0.0_rk

w = (zhigh - zlow)/DBLE(znum-1_ik)

! This is equations (3a) and (3b) from Tauchen (1986)

DO j = 1, znum
	
	z = ZVEC(1) - persistence*ZVEC(j)
	F1 = Normal(z + w/2.0_rk, meaninnov, stdinnov)
	PI(j,1) = F1
	
	DO k = 2, znum - 1
		z = ZVEC(k) - persistence*ZVEC(j)
		F1 = Normal(z + w/2.0_rk, meaninnov, stdinnov)
		F0 = Normal(z - w/2.0_rk, meaninnov, stdinnov)
		PI(j,k) = F1 - F0
	END DO

	z = ZVEC(znum) - persistence*ZVEC(j)
	F0 = Normal(z - w/2.0_rk, meaninnov, stdinnov)
	PI(j,znum) = 1.0_rk - F0

END DO

END SUBROUTINE TAUCHEN2

function normal(x, mean, std)

! program computes the cumulative normal distribution function, using the error function.
!
! the error function (erf.f) is fortran 77 code from http://www.netlib.org/specfun/erf

use kindset

implicit none

real(rk):: x, mean, std, z, normal, derf

intent(in):: x, mean, std

external:: derf

! translate x ~ N(mean, std) into z ~ N(0,1)
z = (x - mean)/std

z = z/dsqrt(2.0_rk)

! compute error function erf(z)
normal = derf(z)

! translate error function into Normal Cumulative Distribution Function
normal = 0.5_rk*(1.0_rk + normal)

end function normal