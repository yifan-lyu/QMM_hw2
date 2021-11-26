PROGRAM InventoriesEdit

! Simplified version of programs to compute the dynamic stochastic equilibrium of the (S,s) inventory model 
! of Khan and Thomas (2007).  While several switches and checks have been removed from these programs for 
! usability, they will reproduce all results in the paper.    
!
! Programs
!
! InventoriesEdit
! Modules
!	KindSet
!	ppslinefit3edit
!	InvToolshededit (InvSstar, InvM, InvPolicy, SetUp, CostFunction, LinSpace, LogSpace)
!	InvWorkShededit (InvInnerE, InvOuter [QfromP, ValueKf, DynEconomyQ])
!
! This set of programs requires the external subroutine DGESV from LAPACK.  
! This is available in the intel math kernel library.  


USE KINDSET
USE ppsplinefit3edit
USE InvToolshededit
USE InvWorkShededit

IMPLICIT NONE

INTEGER(ik):: NKM, NS, NS0, NSM, Nz, Jmax, rknotsV(3), maxVknots, what, numberviolationsKM, &
			  numberviolationsSM, numberviolationsS, numberviolationsQ, numbertry, &
			  Sindex, KMindex, SMindex, md, clmatV, simlength, t, s1, s2, &
			  zi, si, ki, oli, olgo, olgr, ols1, ols2, olsize, olplace, jz, jKM, jSM, jS, &
			  iz1, iz2, ziq
INTEGER(ik), ALLOCATABLE:: numelemV(:,:), ZIsim(:), samebins(:)
REAL(rk):: alpha, thetam, thetan, Atech, delta, abeta, bbeta, masspoint, zibound, &
		   outputcost, Apref, beta, storage, precisionv, precisionm, precisionS, &
		   KMbound(2), Sbound(2), SMbound(2), rhoz, stde, SVECV(2,1), SVECE1(2,1), precisionp, &
		   precisionq, pbounds(2), qbounds(2), distanceKinv, qmiss, rlambda, &
		   distanceInnerV, logbound
REAL(rk), ALLOCATABLE:: Zvec(:,:), PI(:,:), knotsKM(:), knotsS(:), knotsSM(:), hfraction(:), &
						knotsV(:,:), LV(:,:), UV(:,:), dtauV(:), Alphad(:,:), AdjCost(:), &
						csS(:,:), LMATV(:,:), UMATV(:,:,:), DTAUMATV(:,:,:), cV(:,:,:), &
						BetaKM(:,:), BetaSM(:,:), BetaP(:,:), NKagg(:), NMagg(:), &
						BetaPQ(:,:), MUsim(:,:), Ssim(:,:), Ksim(:), Iagg(:), Yagg(:), SB0(:), &
						V(:,:,:), KB(:), QB(:), PB(:), Yd(:,:), MCd(:,:), NCd(:,:), &
						Psim(:), Qsim(:), XR(:,:), SMsim(:), statsPQ(:,:), statsP(:,:), &
						statsK(:,:), statsS(:,:), Sstare(:,:,:), Cagg(:), Xagg(:), &
						Qerror(:), Perror(:), BetaPQInter(:,:), BetaPInter(:,:), &
						BetaKMInter(:,:), BetaSMInter(:,:)

REAL(rk), ALLOCATABLE:: Sprodndist(:,:), hSim(:), SstarSim(:,:)

CHARACTER(30):: datafile, resultfile, indicator, shockfile, initialfile
CHARACTER(30):: OutLoop, RegLoop, OutLoops, RegLoops, InLoops, OutRunLoops, olname


! ********************* !
!!!!!!!!!!!!!!!!!!!!!!!!!
!						!
!	PARAMETER VALUES	!
!						!
!!!!!!!!!!!!!!!!!!!!!!!!!

! Retrieve parameterisation of the model
CALL SETUP(alpha, thetam, thetan, Atech, delta, abeta, bbeta, masspoint, zibound,&
		   outputcost, Apref, beta, storage, precisionv, precisionm, precisionS, precisionp, &
		   precisionq, NKM, KMbound, NS, NS0, Sbound, NSM, SMbound, pbounds, &
		   qbounds, simlength, Jmax, numbertry, datafile, resultfile, shockfile, initialfile)  

! erase the last 4 digits which should be .txt.
olgo = LEN_TRIM(datafile); datafile = datafile(1:olgo-4)
OutLoop(1:3) = 'Out'
RegLoop(1:3) = 'Reg'
olgo = LEN_TRIM(datafile)
OutLoop(4:3+olgo) = datafile
RegLoop(4:3+olgo) = datafile
OutLoop = OutLoop(1:olgo+3)
RegLoop = RegLoop(1:olgo+3)

olgo = LEN_TRIM(OutLoop)
olgr = LEN_TRIM(RegLoop)

! Allocate dimensions
ALLOCATE(Ksim(simlength+1), Psim(simlength), Qsim(simlength), SMsim(simlength+1), &
		 Yd(simlength, Jmax+1), MCd(simlength, Jmax+1), NCd(simlength, Jmax+1), &
		 Cagg(simlength), Xagg(simlength), Iagg(simlength), NKagg(simlength), &
		 NMagg(simlength), AdjCost(simlength), XR(simlength,2), hfraction(simlength), &
		 Yagg(simlength), ZIsim(simlength), QB(simlength), &
		 PB(simlength), KB(simlength+1), SB0(simlength+1))

ALLOCATE(Sprodndist(simlength,Jmax+1), hSim(simlength), SstarSim(simlength,2))
ALLOCATE(Ssim(simlength+1, Jmax+1), MUsim(simlength+1, Jmax+1), Alphad(simlength, Jmax+1), samebins(Jmax))


! ************************************************************************************* !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!																						!
!	Initial Regression Setup: from the baseline or a previous outer loop simulation		!
!																						!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ************************************************************************************* !

! Read data from baseline or previous outer loop simulation

WRITE(*,'(1X,A)', ADVANCE = 'NO') ' Read Previous Outer Simulation? (y = 0) '
READ(*,*) ziq

PRINT*, ziq

IF (ziq.EQ.0) THEN 

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!														!	
	!	 Beginning from a previous Outer Loop Simulation	!
	!														!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	PRINT*, ' Reading previous OutRun '
	
	OPEN(UNIT = 30, FILE = 'OutRun', ACTION = 'READ', STATUS = 'OLD')
	
	READ(30,*) s1

	DO t = 1, simlength
		READ(30,*) ZIsim(t), QB(t), PB(t)
	END DO

	DO t = 1, simlength + 1
		READ(30,*) KB(t), SB0(t)
	END DO
	
	READ(30,*) s2							            ! The iteration number of the saved outerloop

	READ(30, *, END = 100) distanceKinv
	READ(30, *, END = 100) distanceKinv
	READ(30, *, END = 100) distanceKinv
	READ(30,*) Nz
	
	! Discretize the stochastic exogenous productivity process
    ALLOCATE(Zvec(2,Nz), PI(Nz,Nz), statsPQ(Nz,5), statsP(Nz,5),statsK(Nz,5),  statsS(Nz,5), &
             BetaPQ(Nz,3), BetaP(Nz,3), BetaKM(Nz,3), BetaSM(Nz,3), Sstare(NKM, NSM, Nz), &
             BetaPQInter(Nz,3), BetaPInter(Nz,3), BetaKMInter(Nz,3), BetaSMInter(Nz,3), &
             V(NS,NKM*NSM, Nz))

    DO iz1 = 1, Nz
        READ(30,*) Zvec(1,iz1), Zvec(2,iz1)
    END DO
		
    DO iz1 = 1, Nz
        DO iz2 = 1, Nz
            READ(30,*) PI(iz1, iz2)
        END DO
    END DO

	100 CLOSE(30)

	XR(:,1) = DLOG(KB(1:simlength)) 
	XR(:,2) = DLOG(SB0(1:simlength))
	
	WRITE(*,*) ' '	
	
		BetaPQ = 0.0_rk	
		CALL InvPolicy(simlength, 2, Nz, zIsim, DLOG(PB*QB), XR(:,1:2), BetaPQ(1:Nz,1:3), statsPQ)
		WRITE(*,*) ' ' 
		WRITE(*,*) ' PQ Regression Results with SM as a regressor' 
		DO zi = 1, Nz	
			WRITE(*, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaPQ(zi,:)
			WRITE(*,*) '    S.E.      R2      minerror  maxerror     obs.   '
			WRITE(*, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsPQ(zi,:)
		END DO

		BetaP = 0.0_rk
		CALL InvPolicy(simlength, 2, Nz, zIsim, DLOG(PB), XR(:,1:2), BetaP(1:Nz,1:3), statsP)
		WRITE(*,*) ' ' 
		WRITE(*,*) ' P  Regression Results with SM as a regressor ' 
		WRITE(*,*) ' '
		DO zi = 1, Nz	
			WRITE(*, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaP(zi,:)
			WRITE(*,*) '    S.E.      R2      minerror  maxerror     obs.   '
			WRITE(*, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsP(zi,:)
		END DO

		BetaKM = 0.0_rk
		CALL InvPolicy(simlength, 2, Nz, zIsim, DLOG(KB(2:simlength+1)), XR(:,1:2), BetaKM(1:Nz,1:3), statsK)
		WRITE(*,*) ' ' 
		WRITE(*,*) ' KM Regression Results with SM as a regressor ' 
		WRITE(*,*) ' '

		DO zi = 1, Nz	
			WRITE(*, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaKM(zi,:)
			WRITE(*,*) '    S.E.      R2      minerror  maxerror     obs.   '
			WRITE(*, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsK(zi,:)
		END DO

		BetaSM = 0.0_rk
		CALL InvPolicy(simlength, 2, Nz, zIsim, DLOG(SB0(2:simlength+1)), XR(:,1:2), BetaSM(1:Nz,1:3), statsS)
		WRITE(*,*) ' '
		WRITE(*,*) ' SM Regression Results with KM as a regressor ' 
	
		WRITE(*,*) ' '

		DO zi = 1, Nz	
			WRITE(*, '(1X, A)', ADVANCE = 'NO') ' Coefficients  '
			DO t = 1,3
				WRITE(*, '(1X, F10.5)', ADVANCE = 'NO') BetaSM(zi,t)
			END DO
			WRITE(*,*) ' '
			WRITE(*,*) '    S.E.      R2      minerror  maxerror     obs.   '
			WRITE(*, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsS(zi,:)
	    END DO

ELSE

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	!									!
	!	 Beginning from the Baseline	!
	!									!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	PRINT*, ' Reading from baseline ' 
	
	s2 = 0_ik

	OPEN(UNIT = 30, FILE = shockfile, ACTION = 'READ', STATUS = 'OLD')
	
	READ(30,*) t

	IF(t.NE.simlength) THEN
		PRINT*, ' simlength does not match baseline history '
	END IF

	DO t = 1, simlength
		READ(30,*) ZIsim(t), QB(t), PB(t)
	END DO
	
	DO t = 1, simlength + 1
		READ(30,*) KB(t)
	END DO

	! End of file read error results in close of file
	READ(30,*) distanceKinv
	READ(30,*) distanceKinv
	READ(30,*) distanceKinv
	READ(30,*) Nz
	
	! Discretize the stochastic exogenous productivity process
    ALLOCATE(Zvec(2,Nz), PI(Nz,Nz), statsPQ(Nz,5), statsP(Nz,5),statsK(Nz,5),  statsS(Nz,5), &
             BetaPQ(Nz,3), BetaP(Nz,3), BetaKM(Nz,3), BetaSM(Nz,3), Sstare(NKM, NSM, Nz), &
             BetaPQInter(Nz,3), BetaPInter(Nz,3), BetaKMInter(Nz,3), BetaSMInter(Nz,3), &
             V(NS,NKM*NSM, Nz))


    DO iz1 = 1, Nz
        READ(30,*) Zvec(1,iz1), Zvec(2,iz1)
    END DO
		
    DO iz1 = 1, Nz
        DO iz2 = 1, Nz
            READ(30,*) PI(iz1, iz2)
        END DO
    END DO
	
	200 CLOSE(30)
	

	! All regressions are in logs.
	BetaPQ = 0.0_rk
	CALL InvPolicy(simlength, 1, Nz, ZIsim, DLOG(PB*QB), DLOG(KB(1:simlength)), BetaPQ(1:Nz,1:2), statsPQ)

	WRITE(*,*) ' ' 
	WRITE(*,*) ' PQ Initial Regression Results ' 
	WRITE(*,*) ' '

	DO zi = 1, Nz
	
		WRITE(*, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaPQ(zi,:)
		WRITE(*,*) '    S.E.      R2      minerror  maxerror     obs.   '
		WRITE(*, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsPQ(zi,:)
		WRITE(*,*) ' '

	END DO

	BetaP = 0.0_rk
	CALL InvPolicy(simlength, 1, Nz, ZIsim, DLOG(PB), DLOG(KB(1:simlength)), BetaP(1:Nz,1:2), statsP)

	WRITE(*,*) ' ' 
	WRITE(*,*) ' P  Initial Regression Results ' 
	WRITE(*,*) ' '

	DO zi = 1, Nz
	
		WRITE(*, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaP(zi,:)
		WRITE(*,*) '    S.E.      R2      minerror  maxerror     obs.   '
		WRITE(*, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsP(zi,:)
		WRITE(*,*) ' '

	END DO

	BetaKM = 0.0_rk
	CALL InvPolicy(simlength, 1, Nz, ZIsim, DLOG(KB(2:simlength+1)), DLOG(KB(1:simlength)), BetaKM(1:Nz,1:2), statsK)

	WRITE(*,*) ' ' 
	WRITE(*,*) ' KM Initial Regression Results ' 
	WRITE(*,*) ' '

	DO zi = 1, Nz
	
		WRITE(*, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaKM(zi,:)
		WRITE(*,*) '    S.E.      R2      minerror  maxerror     obs.   '
		WRITE(*, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsK(zi,:)
		WRITE(*,*) ' '

	END DO

	BetaSM = 0.0_rk								! Set initial regression for SM to midpoint of bounds
	BetaSM(:,1) = DLOG(SUM(SMbound)/2.0_rk)

END IF


! *********************************************	!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!												!
!   Set Up Period 0 using Steady State Solution !
!												!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ********************************************* !

V = 0.0_rk; distanceInnerV = 0.0_rk

! Get initial conditions for S and MU vector as well as K using the steady state 
! recorded in initialfile.
Ssim = 0.0_rk; MUsim = 0.0_rk; Ksim = 0.0_rk

OPEN(UNIT = 30, FILE = initialfile, ACTION = 'READ', STATUS = 'OLD')	
READ(30,*) s1
DO t = 1, s1
READ(30,*) Ssim(1,t)
END DO

DO t = 1, s1
READ(30,*) MUsim(1,t)
END DO

READ(30,*) Ksim(1)

READ(30,*) indicator

DO si = 1, 2
READ(30,*) SVECV(si,1)
END DO

DO si = 1,2
READ(30,*) SVECE1(si,1)
END DO

READ(30,*, IOSTAT = what, ERR = 300) logbound

CLOSE(30)

! The subroutine comes to line 300 if the file ended after logbound was not available.
300	IF (what.NE.0_ik) THEN 
			WRITE(*,'(1X,A,I4)', ADVANCE = 'YES') ' At the end of the steady state input file, I got an I/O error = ', what
			WRITE(*,*) ' Since logspacing parameter is not available, I will linspace. '
			CLOSE (30)
			WRITE(*,*) ' '
			logbound = 0.0_rk
	END IF

si = LEN_TRIM(indicator)

WRITE(*,'(1X,A, F8.4)', ADVANCE = 'YES') ' Current value of logbound is ', logbound


! ********************************* !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									!
! Grid generation including knots	!
!   SPLINE SETUP FOR 2 FUNCTIONS	!
!									!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ********************************* !

! Generate knots on K, S, SM.  knotsV stacks (S,KM,SM) for each z.
Sindex = 1; KMindex = 2; SMindex = 3

rknotsV = (/NS, NKM, NSM/)

maxVknots = MAXVAL(rknotsV)

ALLOCATE(knotsKM(NKM), knotsS(NS), knotsSM(NSM), knotsV(maxVknots,3))

CALL LINSPACE(KMbound(1), KMbound(2), NKM, knotsKM)
CALL LINSPACE(Sbound(1), Sbound(2), NS, knotsS)
CALL LINSPACE(SMbound(1), SMbound(2), NSM, knotsSM)

! Stick in multivariate spline knots
knotsV = 0.0_rk
knotsV(1:rknotsV(Sindex), Sindex) = knotsS
knotsV(1:rknotsV(KMindex), KMindex) = knotsKM
knotsV(1:rknotsV(SMindex), SMindex) = knotsSM

! Set up univariate splines (for S1)
DEALLOCATE(knotsS, knotsKM, knotsSM); ALLOCATE(knotsS(NS0))
CALL LINSPACE(Sbound(1), Sbound(2), NS0, knotsS)

IF (logbound.GE.precisionv) THEN
	CALL LOGSPACE(knotsS(2)/logbound, Sbound(2), NS0 - 1_ik, knotsS(2:NS0))
	knotsV(1:rknotsV(Sindex), Sindex) = knotsS
END IF


! UNIVARIATE KNOT SETUP.  Most of this follows SplineDriver.F90, in particular the dimension of 
! LV, UV, dtauV and csS should be robust for any one-dimensional example.  
ALLOCATE(LV(1, NS0-3), UV(2, NS0-2), dtauV(NS0-1))
CALL SPLHS(knotsS, NS0-2, indicator, LV, UV, dtauV)
ALLOCATE(csS(4, NS0-1)); csS = 0.0_rk

! Set up multivariate spline (on the space: (S;K,mu))
md = 3_ik

! MULTIVARIATE KNOT SETUP
ALLOCATE(numelemV(4,md))

! preliminary spline fit step, SPFitA computes knot dependent terms, numelem is always (4 x m).
CALL SPFitA0(md, rknotsV, clmatV, numelemV)

! The dimension declared below is robust to general values of m.  However, the one application specific 
! entry is that cV has an additional dimension here, with Nz, because we will interpolate a separate value 
! function for each value of z.  We are replacing a 4-dimensional spline in (S,z,KM,SM) with Nz 
! 3-dimensional splines in (S,KM,SM).  
ALLOCATE(LMATV(md, clmatV), UMATV(2, clmatV+1, md), DTAUMATV(clmatV+1, numelemV(1,md), md), &
		 cV(4, numelemV(3,md), Nz))
cV = 0.0_rk

! first spline fit step: determine function independent LHS terms for slope solution
CALL SPFitA(md, rknotsV, clmatV, knotsV, numelemV, LMATV, UMATV, DTAUMATV)

!!	**********************************	!!
!!  **********************************	!!
!!										!!
!!										!!
!!	ITERATION ON INNER AND OUTER LOOPS	!!		
!!										!!
!!										!!
!!	**********************************	!!
!!  **********************************	!!

InLoops(1:11) = 'InnerSpline'; OutRunLoops(1:6) = 'OutRun'
distanceKinv = 100.0_rk; s1 = s2; rlambda = 0.0_rk		
cV = 0.0_rk; V = 0.0_rk

WRITE(*,'(1X,A,A,A)', ADVANCE = 'NO') ' endpoint condition: ', indicator(1:si), &
    ' skip the Inner Loop? (yes = 0, initial condition = 1) '
READ(*,*) zi

! Read an existing set of spline coefficients either to skip the inner loop or have an initial condition.
IF (zi.EQ.1_ik) THEN 

	WRITE(*,*) ' Reading an initial value function ' 

	OPEN(UNIT = 45, FILE = 'InitialValue', ACTION = 'READ', STATUS = 'OLD')

		DO jz = 1, Nz; DO jKM = 1, NKM; DO jSM = 1, NSM; DO jS = 1, NS
				READ(45,*) V(jS, (jKM -1)*NSM + jSM, jz)
		END DO; END DO; END DO; END DO

	CLOSE (45)

END IF

! Get back to things
ALLOCATE(Qerror(simlength), Perror(simlength))

DO 

	IF (distanceKinv.LE.precisionv*0.75_rk.OR.s1.GT.250_ik) THEN	! Convergence achieved we exit
		EXIT
	ELSE
		s1 = s1 + 1_ik
	END IF
	
	! Create filenames for Outerloop and Regression results which append the iteration number.
	OutLoops(1:olgo) = OutLoop; RegLoops(1:olgr) = RegLoop	
	IF (s1.GT.99) THEN; olsize = 3; ELSEIF (s1.GT.9) THEN; olsize = 2 
	ELSEIF (s1.GT.-1) THEN; olsize = 1; ELSE; EXIT; END IF
	ols1 = s1
	DO oli = 1, olsize
		olplace = 10**(olsize - oli); ols2 = ols1/olplace; olname = ACHAR(ols2+48)	
		OutLoops(olgo+oli:olgo+oli) = olname(1:1); RegLoops(olgr+oli:olgr+oli) = olname(1:1)
		InLoops(11+oli:11+oli) = olname(1:1); OutRunLoops(6+oli:6+oli) = olname(1:1)
		ols1 = ols1 - ols2*olplace	
	END DO

		
		!!!!!!!!!!!!!!!!!!!!!
		!					!		
		!	THE INNER LOOP	!
		!					!
		!!!!!!!!!!!!!!!!!!!!!

IF (zi.NE.0_ik) THEN	! Solve the Inner Loop: call a subroutine


	CALL InvInnerE(abeta, bbeta, masspoint, zibound, outputcost, thetam, thetan, &
				   Atech, beta, Apref, storage, Nz, NS, NS0, NKM, NSM, Zvec, PI,  &
				   knotsS, rknotsV, knotsV, md, numelemV, cV, clmatV, LMATV, UMATV, DTAUMATV, &
				   LV, UV, dtauV, csS, SVECV, SVECE1, indicator, BetaKM, BetaSM, &
				   BetaP, BetaPQ, precisionm, precisionS, precisionv, KMindex, &
				   SMindex, Sindex, V, Sstare, distanceInnerV)

	
	! Catch the splines we've solved, cV and cV(4, numelemV(3,md), Nz)
	OPEN(UNIT = 44, FILE = InLoops, ACTION = 'WRITE', STATUS = 'REPLACE')
	DO zi = 1,Nz; DO ki = 1,numelemV(3,md); DO si = 1,4
		WRITE(44,*) cV(si,ki,zi)
	END DO; END DO; END DO
	WRITE(44,*) logbound
	CLOSE (44)

ELSE	! we have skipped the inner loop, read a stored spline

	OPEN(UNIT = 44, FILE = 'InnerSpline', ACTION = 'READ', STATUS = 'OLD')
	DO zi = 1,Nz; DO ki = 1,numelemV(3,md); DO si = 1,4
		READ(44,*) cV(si,ki,zi)
	END DO; END DO; END DO
	CLOSE (44)

END IF

		!!!!!!!!!!!!!!!!!!!!!
		!					!
		!	THE OUTER LOOP	!
		!					!
		!!!!!!!!!!!!!!!!!!!!!

! Solve the Outer Loop: call a subroutine

WRITE(*,*) ' '

CALL InvOuter(alpha, delta, abeta, bbeta, masspoint, zibound, outputcost, thetam, thetan, &
			  Atech, beta, Apref, storage, Nz, NS, NS0, Zvec, PI, knotsS, &
			  rknotsV, knotsV, md, numelemV, cV, clmatV, LV, UV, dtauV, csS, SVECV, SVECE1, &
			  indicator, BetaKM, BetaSM,  BetaP, BetaPQ, &
			  precisionm, precisionS, precisionp, precisionq, simlength, zIsim, pbounds, &
			  qbounds, KMbound, SMbound, Jmax, 20_ik, s1, PB, MUsim, Ssim, Ksim, Qsim, &
			  Psim, qmiss, Yd, MCd, NCd, Cagg, Yagg, Xagg, Iagg, NKagg, NMagg, Alphad, &
			  AdjCost, distanceKinv, Qerror, Perror, numberviolationsKM, numberviolationsSM, &
			  numberviolationsS, numberviolationsQ, Sprodndist)



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!							!
		!	 THE REGRESSION STEP	!
		!							!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SMsim = SUM(Ssim*MUsim, DIM = 2) 
XR(:,1) = DLOG(Ksim(1:simlength))
XR(:,2) = DLOG(SMsim(1:simlength))

	BetaPQInter = 0.0_rk
	CALL InvPolicy(simlength, 2, Nz, zIsim, DLOG(Psim*Qsim), XR(:,1:2), BetaPQInter(1:Nz,1:3), statsPQ)
	WRITE(*,*) ' ' 
	WRITE(*,*) ' PQ Regression Results with SM as regressor ' 
	WRITE(*,*) ' '
	DO zi = 1, Nz	
		WRITE(*, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaPQInter(zi,1:3)
		WRITE(*,*) '    S.E.      R2      minerror  maxerror     obs.   '
		WRITE(*, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsPQ(zi,:)
	END DO	

	
	BetaPQ = rlambda*BetaPQ + (1.0_rk - rlambda)*BetaPQInter(1:Nz,1:3)
	
	BetaPInter = 0.0_rk
	CALL InvPolicy(simlength, 2, Nz, zIsim, DLOG(Psim), XR(:,1:2), BetaPInter(1:Nz,1:3), statsP)
	WRITE(*,*) ' ' 
	WRITE(*,*) ' P  Regression Results ' 
	WRITE(*,*) ' '
	DO zi = 1, Nz	
		WRITE(*, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaPInter(zi,1:3)
		WRITE(*,*) '    S.E.      R2      minerror  maxerror     obs.   '
		WRITE(*, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsP(zi,:)
	END DO
	
	BetaP = rlambda*BetaP + (1.0_rk - rlambda)*BetaPInter(1:Nz,1:3)	

	BetaKMInter = 0.0_rk
	CALL InvPolicy(simlength, 2, Nz, zIsim, DLOG(Ksim(2:simlength+1)), XR(:,1:2), BetaKMInter(1:Nz,1:3), statsK)
	WRITE(*,*) ' ' 
	WRITE(*,*) ' KM Regression Results ' 
	WRITE(*,*) ' '

	DO zi = 1, Nz	
		WRITE(*, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaKMInter(zi,1:3)
		WRITE(*,*) '    S.E.      R2      minerror  maxerror     obs.   '
		WRITE(*, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsK(zi,:)
	END DO

	
	BetaKM = rlambda*BetaKM +  (1.0_rk - rlambda)*BetaKMInter(1:Nz,1:3)

	BetaSMInter = 0.0_rk
	
		CALL InvPolicy(simlength, 2, Nz, zIsim, DLOG(SMsim(2:simlength+1)), XR(:,1:2), BetaSMInter(1:Nz,1:3), statsS)
		WRITE(*,*) ' ' 
		WRITE(*,*) ' SM Regression Results ' 


	WRITE(*,*) ' '

	DO zi = 1, Nz	
		WRITE(*, '(1X, A)', ADVANCE = 'NO') ' Coefficients  '
		DO t = 1,3
			WRITE(*, '(1X, F10.5)', ADVANCE = 'NO') BetaSMInter(zi,t)
		END DO
		WRITE(*,*) ' '
		WRITE(*,*) '    S.E.      R2      minerror  maxerror     obs.   '
		WRITE(*, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsS(zi,:)
	END DO

	IF (BetaSM(1,2).NE.0.0_rk) THEN	! handles first iterations using baseline data where BetaSM(:,2:3) = 0.0
		BetaSM = rlambda*BetaSM + (1.0_rk - rlambda)*BetaSMInter(:,1:3)
	ELSE 
		BetaSM = BetaSMInter(:,1:3)
	END IF 

distanceKinv = MAXVAL(DABS(Ksim - KB))

KB = Ksim
PB = Psim
QB = Qsim
SB0 = SMsim

WRITE(*,*) ' ' 

WRITE(*,FMT = '(1X,A,I3,2(A,F8.4))') ' overall iteration ', s1, ' dInnverV: ', distanceInnerV, ' dCapital: ', distanceKinv

WRITE(*,'(1X,A,F8.4, F8.4)',ADVANCE = 'YES') ' Maximum Q and P errors: ', MAXVAL(DABS(Qerror)), MAXVAL(DABS(Perror))
WRITE(*,*) ' ' 

WRITE(*,'(1X,A,F8.4, F8.4)',ADVANCE = 'YES') ' Average Q and P errors: ', SUM(DABS(Qerror))/simlength, SUM(DABS(Perror))/simlength
WRITE(*,*) ' ' 

WRITE(*,'(1X,A,I5)', ADVANCE = 'YES') 'KM grid violations in simulation: ', numberviolationsKM
WRITE(*,'(1X,A,I5)', ADVANCE = 'YES') 'SM grid violations in simulation: ', numberviolationsSM
WRITE(*,'(1X,A,I5)', ADVANCE = 'YES') ' S grid violations in simulation: ', numberviolationsS
WRITE(*,'(1X,A,I5)', ADVANCE = 'YES') ' Q grid violations in simulation: ', numberviolationsQ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									!		
!		SAVE RUNNING RESULTS		!
!									!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Save Results

OPEN(UNIT = 45, FILE = OutLoops, ACTION = 'WRITE', STATUS = 'REPLACE')	

WRITE(45,*) datafile
WRITE(45,*) s1
WRITE(45,*) simlength
WRITE(45,*) Jmax
WRITE(45,*) outputcost
DO t = 1, simlength+1; WRITE(45,*) Ksim(t); END DO
DO t = 1, simlength; WRITE(45,*) Psim(t); END DO
DO t = 1, simlength; WRITE(45,*) Qsim(t); END DO
DO t = 1, simlength; WRITE(45,*) Cagg(t); END DO
DO t = 1, simlength; WRITE(45,*) Yagg(t); END DO
DO t = 1, simlength; WRITE(45,*) Xagg(t); END DO
DO t = 1, simlength; WRITE(45,*) Iagg(t); END DO
DO t = 1, simlength; WRITE(45,*) NKagg(t); END DO

DO t = 1, simlength+1; DO zi = 1, Jmax; WRITE(45,*) Ssim(t,zi); END DO; END DO
DO t = 1, simlength+1; DO zi = 1, Jmax; WRITE(45,*) MUsim(t,zi); END DO; END DO

DO t = 1, simlength; DO zi = 1, Jmax + 1; WRITE(45,*) Yd(t,zi); END DO; END DO
DO t = 1, simlength; DO zi = 1, Jmax + 1; WRITE(45,*) MCd(t,zi); END DO; END DO
DO t = 1, simlength; DO zi = 1, Jmax + 1; WRITE(45,*) NCd(t,zi); END DO; END DO

DO t = 1, simlength; DO zi = 1, Jmax; WRITE(45,*) Alphad(t,zi); END DO; END DO

DO t = 1, simlength; WRITE(45,*) AdjCost(t); END DO

DO t = 1, simlength; WRITE(45,*) hfraction(t); END DO

! Parameters
WRITE(45,*) ' '
WRITE(45,*) ' '

WRITE(45,*) alpha; WRITE(45,*) thetam; WRITE(45,*) thetan
WRITE(45,*) Atech; WRITE(45,*) rhoz 
WRITE(45,*) stde; WRITE(45,*) Nz; WRITE(45,*) abeta; WRITE(45,*) bbeta
WRITE(45,*) masspoint; WRITE(45,*) zibound; WRITE(45,*) outputcost; WRITE(45,*) beta
WRITE(45,*) Apref; WRITE(45,*) (Sbound(t), t = 1,2); WRITE(45,*) NS0, NS
WRITE(45,*) (KMbound(t),t = 1,2); WRITE(45,*) NKM
WRITE(45,*)(SMbound(t), t = 1,2); WRITE(45,*) NSM
WRITE(45,*) precisionv; WRITE(45,*) precisionm; WRITE(45,*) precisionS 
WRITE(45,*) precisionp; WRITE(45,*) precisionq; WRITE(45,*) simlength
WRITE(45,*) Jmax; WRITE(45,*) (pbounds(t), t = 1,2)
WRITE(45,*) (qbounds(t), t = 1,2); WRITE(45,*) resultfile
WRITE(45,*) shockfile; WRITE(45,*) initialfile; WRITE(45,*) indicator
WRITE(45,*) (SVECV(t,1), t = 1,2); WRITE(45,*) (SVECE1(t,1),t = 1,2)

WRITE(45,*) ' '; WRITE(45,*) ' '

WRITE(45,*) ' ' 
WRITE(45,*) ' ITERATION ', s1, ' with changes ', distanceKinv
WRITE(45,*) ' '

WRITE(45,*) ' PQ Regression Results with SM as regressor '; DO zi = 1, Nz	
WRITE(45, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaPQInter(zi,:)
WRITE(45,*) '    S.E.      R2      minerror  maxerror     obs.   '
WRITE(45, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsPQ(zi,:); END DO

WRITE(45,*) ' ' 
WRITE(45,*) ' P  Regression Results '; DO zi = 1, Nz	
WRITE(45, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaPInter(zi,:)
WRITE(45,*) '    S.E.      R2      minerror  maxerror     obs.   '
WRITE(45, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsP(zi,:); END DO
WRITE(45,*) ' ' 
WRITE(45,*) ' KM Regression Results '; DO zi = 1, Nz	
WRITE(45, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaKMInter(zi,:)
WRITE(45,*) '    S.E.      R2      minerror  maxerror     obs.   '
WRITE(45, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsK(zi,:); END DO
WRITE(45,*) ' ' 
WRITE(45,*) ' SM Regression Results '; DO zi = 1, Nz	
WRITE(45, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaSMInter(zi,:)
WRITE(45,*) '    S.E.      R2      minerror  maxerror     obs.   '
WRITE(45, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsS(zi,:); END DO

CLOSE(45)

! SAVE CURRENT REGRESSION RESULTS 

OPEN(UNIT = 46, FILE = RegLoops, ACTION = 'WRITE', STATUS = 'REPLACE')	
WRITE(46,*) ' ' 

WRITE(46,'(1X, A, A, I4, 2(A, F10.6))', ADVANCE = 'YES') datafile, ' ITERATION ', s1, ' dInnerV ', &
        distanceInnerV, ' with changes in K ', distanceKinv
WRITE(46,*) ' ' 

WRITE(46,'(1X,A,F8.6)') ' rlambda weight on existing regressions ', rlambda
WRITE(46,*) ' '
WRITE(46,*) ' ' 

WRITE(46,*) ' PQ Regression Results '; DO zi = 1, Nz	
WRITE(46, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaPQInter(zi,:)
WRITE(46,*) '    S.E.      R2      minerror  maxerror     obs.   '
WRITE(46, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsPQ(zi,:); END DO

WRITE(46,*) ' ' 
WRITE(46,*) ' P  Regression Results '; DO zi = 1, Nz	
WRITE(46, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaPInter(zi,:)
WRITE(46,*) '    S.E.      R2      minerror  maxerror     obs.   '
WRITE(46, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsP(zi,:); END DO
WRITE(46,*) ' ' 
WRITE(46,*) ' KM Regression Results '; DO zi = 1, Nz	
WRITE(46, '(1X, A, 3(F10.5))', ADVANCE = 'YES') ' Coefficients  ', BetaKMInter(zi,:)
WRITE(46,*) '    S.E.      R2      minerror  maxerror     obs.   '
WRITE(46, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsK(zi,:); END DO
WRITE(46,*) ' ' 
WRITE(46,*) ' K*SM or SM Regression Results ' 
WRITE(46,*) ' '

DO zi = 1, Nz	
	WRITE(46, '(1X, A)', ADVANCE = 'NO') ' Coefficients  '
	DO t = 1,3
		WRITE(46, '(1X, F10.5)', ADVANCE = 'NO') BetaSMInter(zi,t)
	END DO
	WRITE(46,*) ' '
	WRITE(46,*) '    S.E.      R2      minerror  maxerror     obs.   '
	WRITE(46, '(1X,4(F10.5), F10.0)', ADVANCE = 'YES') statsS(zi,:)
END DO


WRITE(46,*) ' ' 
WRITE(46,*) ' ' 
WRITE(46,'(1X,A,F8.4, F8.4)', ADVANCE = 'YES') ' Maximum Q and P errors: ', &
MAXVAL(DABS(Qerror)), MAXVAL(DABS(Perror))
WRITE(46,*) ' ' 
WRITE(46,'(1X,A,F8.4, F8.4)', ADVANCE = 'YES') ' Average Q and P errors: ', &
SUM(DABS(Qerror))/simlength, SUM(DABS(Perror))/simlength
WRITE(46,*) ' ' 
WRITE(46,*) ' '
WRITE(46,'(1X,A,I5)', ADVANCE = 'YES') 'KM grid violations in simulation: ', numberviolationsKM
WRITE(46,'(1X,A,I5)', ADVANCE = 'YES') 'SM grid violations in simulation: ', numberviolationsSM
WRITE(46,'(1X,A,I5)', ADVANCE = 'YES') ' S grid violations in simulation: ', numberviolationsS
WRITE(46,'(1X,A,I5)', ADVANCE = 'YES') ' Q grid violations in simulation: ', numberviolationsQ
CLOSE(46)

! SAVE RUNNING RESULTS FOR FUTURE USE

OPEN(UNIT = 47, FILE = OutRunLoops, ACTION = 'WRITE', STATUS = 'REPLACE')	
WRITE(47,*) simlength

DO t = 1, simlength
WRITE(47,*) zIsim(t), QB(t), PB(t)
END DO

DO t = 1, simlength + 1
WRITE(47,*) KB(t), SB0(t)
END DO

WRITE(47,*) s1

WRITE(47,*) Nz          ! 3 dummy write statements for backwards compatibility
WRITE(47,*) Nz
WRITE(47,*) Nz
WRITE(47,*) Nz

DO iz1 = 1, Nz
    WRITE(47,*) Zvec(1,iz1), Zvec(2,iz1)
END DO
		
DO iz1 = 1, Nz
    DO iz2 = 1, Nz
        WRITE(47,*) PI(iz1, iz2)
    END DO
END DO

CLOSE(47)

END DO		! THIS IS THE END OF THE ITERATIVE INNER/OUTER LOOP PART

WRITE(*,'(1X,A)', ADVANCE = 'NO' ) ' The Dynamic Inventory programs using ppsplinefit3 have completed.  Press enter. '
READ(*,*)

END PROGRAM InventoriesEdit