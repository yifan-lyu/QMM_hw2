MODULE InvWorkShededit

! Module with inner and outer loop for the inventory model of Khan and Thomas (2007).

USE KindSet
USE ppsplinefit3edit
USE InvToolShededit

IMPLICIT NONE

PRIVATE

PUBLIC:: InvInnerE, InvOuter, QfromP, ValueKf, DynEconomyQ

CONTAINS

!!	******************************************************************************************	!!
!!																								!!
!!																								!!
!!											InvInner											!!
!!																								!!
!!																								!!
!!	******************************************************************************************	!!


SUBROUTINE InvInnerE(abeta, bbeta, masspoint, zibound, outputcost, thetam, thetan, &
					 Atech, beta, Apref, storage, Nz, NS, NS0, NKM, NSM, Zvec, PI, knotsS, rknotsV, &
					 knotsV, md, numelemV, cV, clmatV, LMATV, UMATV, DTAUMATV, LV, UV, dtauV, csV, &
					 SVECV, SVECE1, indicator, BetaKM, BetaSM, BetaP, BetaPQ, &
					 precisionm, precisionS, precisionv, KMindex, SMindex, Sindex, V, Sstare, &
					 distanceInnerV)

INTEGER(ik):: Nz, NS, NS0, NKM, NSM, md, jz, jS, jKM, jSM, numelemV(4,md), clmatV, jzf, &
			  rknotsV(md), KMindex, SMindex, Sindex, s1, extrapolatek, &
			  extrapolates, accelerate

REAL(rk):: thetam, thetan, Atech, beta, Apref, storage, knotsS(NS0), &
		   knotsV(clmatV+2,md), LMATV(md, clmatV), UMATV(2, clmatV+1, md), costterm, &
		   DTAUMATV(clmatV+1, numelemV(1,md), md), zval, gzval, SMval, KMval, SVECV(2,1), SVECE1(2,1), &
		   cV(4, numelemV(3,md), Nz), LV(1, NS0-3), UV(2, NS0-2), dtauV(NS0-1), KMf, SMf, Phat, &
		   PQhat, cysol, BetaKM(Nz,3), BetaSM(Nz,3), BetaP(Nz,3), BetaPQ(Nz,3), &
		   V(NS, NKM*NSM, Nz), nsol, TV(NS, NKM*NSM, Nz), V1zf(Nz), V1(NS0), alphai, &
		   E1(NS0), zitilde, msol, precisionv, precisionm, precisionS, Zvec(2,Nz),  &
		   PI(Nz, Nz), distanceInnerV, cSV(4,NS0-1), csE1(4, NS0-1), Sval(1), Ea, &
		   abeta, bbeta, masspoint, zibound, outputcost, TV1(NS), Smax, Smin, gooble, gooble2, &
		   E1c(1,1), deltaKM, deltaSM, KMfuse, SMfuse, Vfuse, DVfkuse, DVfsuse, csTV1(4,NS-1), &
		   zitilderaw, Sstar, Sstare(Nz, NKM*NSM), TSstare(Nz, NKM*NSM), &
		   AlphaStare(Nz, NKM*NSM, NS), ZiRawStare(Nz, NKM*NSM, NS), MStare(Nz, NKM*NSM, NS), &
		   TAlphaStare(Nz, NKM*NSM, NS), TZiRawStare(Nz, NKM*NSM, NS), TMStare(Nz, NKM*NSM, NS), &
		   TE1Stare(Nz, NKM*NSM, NS), E1Stare(Nz, NKM*NSM, NS), PQhatStare(Nz,NKM*NSM), &
		   XI, TXI(Nz, NKM*NSM, NS), gza, nsola, cysola, rwsola, Sfa(1), begin, finish

CHARACTER(20):: indicator 

INTENT(IN):: thetam, thetan,Atech, beta, Apref, storage, Nz, NS, NS0, NKM, NSM, &
			 Zvec, PI, knotsS, rknotsV, knotsV, md, numelemV, clmatV, LMATV, UMATV, &
			 DTAUMATV, LV, UV, dtauV, SVECV, SVECE1, indicator, BetaKM, BetaSM, BetaP, BetaPQ, &
			 precisionm, precisionS, precisionv, abeta, bbeta, masspoint, zibound, &
			 outputcost, Sindex, SMindex, KMindex
			 
INTENT(OUT):: Sstare, distanceInnerV
INTENT(INOUT):: cV, csV, V


! The Inner Loop, requires V, initial function, and also BetaKM, BetaSM, BetaP and BetaQ
! policy functions that are linear in (/1.0, KM, SM/) and used to forecast future states
! and prices.
s1 = 0_ik; distanceInnerV = 500.0_rk; Sstare = -2.0_rk
accelerate = 0_ik

CALL CPU_TIME(begin)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!																	!	
!	Initial Value Function assuming no future value of inventories	!
!																	!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (MAXVAL(V).EQ.0.0_rk) THEN

	IF (MAXVAL(DABS(cV)).GT.0.0_rk) THEN		
		s1 = 1_ik			
	END IF

DO jz = 1, Nz								! THE 
	zval = Zvec(1,jz)
	gzval = Zvec(2, jz)

	DO jKM = 1, NKM							! AGGREGATE
		KMval = knotsV(jKM,KMindex)

		DO jSM = 1, NSM						! STATE
			
			SMval = knotsV(jSM,SMindex)

			Phat = DOT_PRODUCT(BetaP(jz,1:3),  (/1.0_rk, DLOG(KMval), SMval/)); Phat = DEXP(Phat)	
	
			DO jS = 1,NS

				Sval(1) = knotsV(jS,Sindex)

				IF (s1.EQ.0_ik) THEN	! We have no spline coefficients we might use toward a V0

					gooble = Phat*gzval*Atech
					gooble2 = gooble*thetan*(Sval(1)**thetam)/Apref
					gooble2 = gooble2**(1.0_rk/(1.0_rk - thetan))

					V(jS, (jKM -1)*NSM + jSM, jz) = gooble*(Sval(1)**thetam)*(gooble2**thetan) - Apref*gooble2

					
				ELSE						! We have spline coefficients we might use toward a V0


					V(jS, (jKM -1)*NSM + jSM, jz) = FastSplinEval(md, numelemV(3,md), cV(:,:,jz), rknotsV-1_ik, clmatV+2_ik, &
																  knotsV, (/Sval, KMval, SMval/), (/0_ik, 0_ik, 0_ik/))

				END IF


			END DO	! S

		END DO		! SM

	END DO			! KM

END DO				! z


END IF

AlphaStare=0.0_rk; ZiRawStare=0.0_rk; MStare=0.0_rk; TAlphaStare=0.0_rk; TZiRawStare=0.0_rk; TMStare = 0.0_rk
Sstare = 0.0_rk; TSstare = 0.0_rk

s1 = 0_ik


DO	! Contraction Mapping on V(S;KM,SM,z)		

	IF (distanceInnerV.LT.precisionv/10.0_rk.OR.s1.GT.500_ik) THEN
		EXIT
	ELSE
		s1 = s1 + 1_ik
				
		    IF (s1.GT.10_ik.AND.accelerate.LT.9_ik) THEN
		        accelerate = accelerate + 1_ik
            ELSE
                accelerate = 0_ik
            END IF		       
        
	END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!										!
!	Update the Multivariate Spline		!
!										!		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This call has been changed for use with ppslinefit3
DO jz = 1, Nz
	CALL SPFitB(md, clmatV, rknotsV, numelemV, LMATV, UMATV, DTAUMATV, V(:,:,jz), cV(:,:, jz))
END DO

DO jz = 1, Nz							! THE 
	zval = Zvec(1,jz)
    gzval = Zvec(2, jz)
    
DO jKM = 1, NKM							! AGGREGATE
	KMval = knotsV(jKM,KMindex)

DO jSM = 1, NSM							! STATE
	SMval = knotsV(jSM,SMindex)
	
	deltaKM = 0.0_rk
	deltaSM = 0.0_rk
	extrapolatek = 0_ik
	extrapolates = 0_ik
	KMfuse = 0.0_rk
	SMfuse = 0.0_rk

	! forecast future state
	
	KMf =  DOT_PRODUCT(BetaKM(jz,1:3), (/1.0_rk, DLOG(KMval), DLOG(SMval)/)); KMf = DEXP(KMf)
	Phat = DOT_PRODUCT(BetaP(jz,1:3),  (/1.0_rk, DLOG(KMval), DLOG(SMval)/)); Phat = DEXP(Phat)
	SMf =  DOT_PRODUCT(BetaSM(jz,1:3), (/1.0_rk, DLOG(KMval), DLOG(SMval)/)); SMf = DEXP(SMf)
	PQhat = DOT_PRODUCT(BetaPQ(jz,1:3),  (/1.0_rk, DLOG(KMval), DLOG(SMval)/)); PQhat = DEXP(PQhat)
	
	IF (KMf.GE.knotsV(NKM,KMindex)) THEN ! upper knots violation for KMf
		WRITE(*,'(1X,A,F12.3, A, 4(F12.3))', ADVANCE = 'YES') ' KMf = ', KMf, '   (z,gz,K,SM) = ', zval, gzval, KMval, SMval
		deltaKM = KMf - knotsV(NKM,KMindex); extrapolatek = 1; KMfuse = knotsV(NKM,KMindex)

	ELSEIF (KMf.LT.knotsV(1,KMindex)) THEN
		WRITE(*,'(1X,A,F12.3, A, 4(F12.3))', ADVANCE = 'YES') ' KMf = ', KMf, '   (z,gz,K,SM) = ', zval, gzval, KMval, SMval
		deltaKM = KMf - knotsV(1,KMindex); extrapolatek = 1; KMfuse = knotsV(1,KMindex)
	END IF

	IF (SMf.GE.knotsV(NSM,SMindex)) THEN  
		WRITE(*,'(1X,A,F12.3, A, 4(F12.3))', ADVANCE = 'YES') ' SMf = ', SMf, '   (z,gz,K,SM) = ', zval, gzval, KMval, SMval
		deltaSM = SMf - knotsV(NSM,SMindex); extrapolates = 1; SMfuse = knotsV(NSM,SMindex)

	ELSEIF (SMf.LT.knotsV(1,SMindex)) THEN
		WRITE(*,'(1X,A,F12.3, A, 4(F12.3))', ADVANCE = 'YES') ' SMf = ', SMf, '   (z,gz,K,SM) = ', zval, gzval, KMval, SMval
		deltaSM = SMf - knotsV(1,SMindex); extrapolates = 1; SMfuse = knotsV(1,SMindex)

	END IF
	
	!!!!!!!!!!!!!!!!!
	!				!
	!	The V STEP	!
	!				!	
	!!!!!!!!!!!!!!!!!

	! Note the 3 step procedure below.  First we collapse V to a univariate spline.  Next we use this 
	! in InvM to retrieve a grid of E1 values.  We fit a univariate spline to this E1.  Third, we 
	! maximize this E1 to determine Ea and S*.

	! costterm is Apref (which is wage*Phat) if adjustment costs are in units of time, if they are in 
	! units of output, it is simply Phat.  
	costterm = (1.0_rk - outputcost)*Apref + outputcost*Phat	

	! fit a univariate spline to V[zi](S;KM,SM) having fixed (KM,SM) and zi (the current z(i) index)

	DO jS = 1,NS0
		Sval = knotsS(jS)
	
		DO jzf = 1,Nz
			
			IF (extrapolatek.EQ.0_ik.AND.extrapolates.EQ.0_ik) THEN ! no knot violations
				V1zf(jzf) = FastSplinEval(md, numelemV(3,md), cV(:,:,jzf), rknotsV-1_ik, clmatV+2_ik, knotsV, &
										  (/Sval, KMf, SMf/), (/0_ik, 0_ik, 0_ik/))

			ELSEIF (extrapolatek.EQ.0_ik.AND.extrapolates.EQ.1_ik) THEN		! only SM violation
					Vfuse = FastSplinEval(md, numelemV(3,md), cV(:,:,jzf), rknotsV-1_ik, clmatV+2_ik, knotsV, &
										  (/Sval, KMf, SMfuse/), (/0_ik, 0_ik, 0_ik/))
					DVfsuse = FastSplinEval(md, numelemV(3,md), cV(:,:,jzf), rknotsV-1_ik, clmatV+2_ik, knotsV, &
										  (/Sval, KMf, SMfuse/), (/0_ik, 0_ik, 1_ik/))
					V1zf(jzf) = Vfuse + deltaSM*DVfsuse

			ELSEIF (extrapolatek.EQ.1_ik.AND.extrapolates.EQ.0_ik) THEN		! only KM violation
					Vfuse = FastSplinEval(md, numelemV(3,md), cV(:,:,jzf), rknotsV-1_ik, clmatV+2_ik, knotsV, &
										  (/Sval, KMfuse, SMf/), (/0_ik, 0_ik, 0_ik/))
					DVfkuse = FastSplinEval(md, numelemV(3,md), cV(:,:,jzf), rknotsV-1_ik, clmatV+2_ik, knotsV, &
										  (/Sval, KMfuse, SMf/), (/0_ik, 1_ik, 0_ik/))
					V1zf(jzf) = Vfuse + deltaKM*DVfkuse

			ELSE															! violation on both SM and KM	
				    Vfuse = FastSplinEval(md, numelemV(3,md), cV(:,:,jzf), rknotsV-1_ik, clmatV+2_ik, knotsV, &
										  (/Sval, KMfuse, SMfuse/), (/0_ik, 0_ik, 0_ik/))
					DVfsuse = FastSplinEval(md, numelemV(3,md), cV(:,:,jzf), rknotsV-1_ik, clmatV+2_ik, knotsV, &
										  (/Sval, KMfuse, SMfuse/), (/0_ik, 0_ik, 1_ik/))
					DVfkuse = FastSplinEval(md, numelemV(3,md), cV(:,:,jzf), rknotsV-1_ik, clmatV+2_ik, knotsV, &
										  (/Sval, KMfuse, SMfuse/), (/0_ik, 1_ik, 0_ik/))
					V1zf(jzf) = Vfuse + deltaSM*DVfsuse + deltaKM*DVfkuse
			END IF
			
			V1zf(jzf) = PI(jz, jzf)*V1zf(jzf)

		END DO
	
		V1(jS) = SUM(V1zf)

	END DO

	! Now we fit a spline to the value of taking out any S level at the end of the period.
	CALL SPpp(V1, NS0 -2_ik, 1, LV, UV, dtauV, indicator, SVECV, csV(:,:))
	! This converts the multivariate V to a univariate V1, csS is now the coefficients
	! to the beginning of next period value function E[V(S';z(j),Kf,SMf)|z(i)].
	
	!!  *****************************************************************************************************  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                                                                           !                                                                                                                                                                                                                                           
	!   Solve the Optimisation problem:                                                                         !
	!                                                                                                           !    
	!       0. Use csV and InvMSUP to obtain msol, stored as TMStare(jz, (jKM-1)*NSM + jSM, jS), and csE1.      !
	!       1. Bracket S* in L New Examination of Objective                                                     !
	!       2. Call InvSstarSup to obtain Sstar and Ea (stored in TSstare(jz, (jKM -1)*NSM + jSM)               !
	!       3. Compute alphai, stored as TAlphaStare(jz, (jKM-1)*NSM + jSM, jS), to update TV1                  !
	!          stored asTV(1:NS, (jKM -1)*NSM + jSM, jz).                                                       !
	!                                                                                                           !    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!  *****************************************************************************************************  !!   
	
	IF (accelerate.EQ.0_ik) THEN
	    
	    ! Solve the optimisation problem
	    
	    DO jS = 1, NS0		! Looping over values for mid-period inventory holdings, S1

		    Sval(1) = knotsS(jS)
		    ! We now construct a function E1(S) using the new univariate function V1.
		    CALL InvM(NS0, thetam, thetan, Atech, Phat, beta, precisionm, Apref, storage, 1.0_rk, knotsS, &
                      csV, (/0.0_rk, Sval(1)/), gzval, Sval(1), msol, nsol, cysol, E1(jS))

		    TMStare(jz, (jKM-1)*NSM + jSM, jS) = msol
		    TE1Stare(jz, (jKM-1)*NSM + jSM, jS) = E1(jS)
	    END DO	! Ends the loop over grid values for S1


	    ! In this 3rd step, we fit a spline, csS to E1 of equation (12), the mid-period value function.   
	    ! We use S*(KM, SM, zval), first fit a spline to E1, csS is now the spline coefficients to E1.
	    CALL SPpp(E1, NS0-2_ik, 1, LV, UV, dtauV, indicator, SVECE1, csE1(:,:))
    	
	    CALL InvSstar(NS0, PQhat, precisionS, knotsS, csE1, (/knotsS(1), knotsS(NS0)/), Sstar, Ea)

	    Sval(1) = Sstar

	    TSstare(jz, (jKM -1)*NSM + jSM) = Sstar
	    PQhatStare(jz, (jKM -1)*NSM + jSM) = PQhat 

	    ! update V
	    alphai = 0.0_rk	

    DO jS = 1, NS

	    ! determine the thresholds using Ea and E1
	    Sval(1) = knotsV(jS,Sindex)

	    CALL SPeval(csE1, knotsS, NS0-2_ik, 1, Sval, 1, 0, E1c) 

	    zitilderaw = (PQhat*Sval(1) - E1c(1,1) + Ea)/costterm
        zitilde = DMIN1(DMAX1(0.0_rk, zitilderaw),zibound)

	    ! Generalized cost function to deliver probability because integer 5th argument is included.
	    alphai = CostFunction(zitilde, abeta, bbeta, masspoint, zibound, 0_ik)  
	    XI = CostFunction(zitilde, abeta, bbeta, masspoint, zibound)
	    TV1(jS) = alphai*(PQhat*Sval(1) + Ea) - costterm*XI

	    TV1(jS) = TV1(jS) + (1.0_rk - alphai)*E1c(1,1)
	    TAlphaStare(jz, (jKM-1)*NSM + jSM, jS) = alphai
	    TXI(jz, (jKM-1)*NSM + jSM, jS) = XI
	    TZiRawStare(jz, (jKM-1)*NSM + jSM, jS) = zitilderaw

    END DO
    
ELSE
    
    ! As accelerate is not equal to 0, we do not solve the optimisation problem, 
    ! instead we use existing decision rules to update the value function
    
    gza = gzval		                        ! c production function construction
    gza = Atech*gza		
    
    DO jS = 1, NS
    
        Sval(1) = knotsS(jS)
        msol = TMStare(jz, (jKM-1)*NSM + jSM, jS)        
        Sfa(1) = Sval(1) - msol
        CALL SPeval(csV, knotsS, NS0 - 2_ik, 1, Sfa, 1, 0, E1c)
        nsola = thetan*Phat*gza*(msol**thetam)/Apref
        nsola = nsola**(1.0_rk/(1.0_rk - thetan))
        cysola = gza*(msol**thetam)*(nsola**thetan)  - storage*Sfa(1)
        rwsola = Phat*cysola - Apref*nsola
        E1(jS) = rwsola + beta*E1c(1,1)
        
    END DO        
    
    CALL SPpp(E1, NS0-2_ik, 1, LV, UV, dtauV, indicator, SVECE1, csE1(:,:))
    
    Sstar = TSstare(jz, (jKM -1)*NSM + jSM)
    Sval(1) = Sstar
    CALL SPeval(csE1, knotsS, NS0 - 2_ik, 1, Sval, 1, 0, E1c)
    Ea = -PQhat*Sval(1) + E1c(1,1)
        
    DO jS = 1, NS

	    Sval(1) = knotsV(jS,Sindex)	    
	    
	    CALL SPeval(csE1, knotsS, NS0-2_ik, 1, Sval, 1, 0, E1c) 	    
	    
	    alphai = TAlphaStare(jz, (jKM-1)*NSM + jSM, jS)
	    XI = TXI(jz, (jKM-1)*NSM + jSM, jS)

	    TV1(jS) = alphai*(PQhat*Sval(1) + Ea) - costterm*XI
	    TV1(jS) = TV1(jS) + (1.0_rk - alphai)*E1c(1,1)	   

    END DO
    
END IF          ! The conditional decision to optimize or skip the computation of decision rules    
    
! See TW above for more information on why we store the spline data in this manner.
TV(1:NS, (jKM -1)*NSM + jSM, jz) = TV1

! NS MUST EQUAL NS0 for this to be sensible.
CALL SPpp(TV1, NS -2_ik, 1, LV, UV, dtauV, indicator, SVECE1, csTV1(:,:))


END DO				! THE 

END DO				! STATE

END DO				! STOPS HERE

distanceInnerV = MAXVAL(DABS(TV-V))
Smin = MINVAL(TSstare); Smax = MAXVAL(TSstare)

WRITE(*,*) ' '

IF (accelerate.EQ.0_ik) THEN

WRITE(*, FMT = '(1X, A, I3, A, F6.2, A, F6.2, A, E9.3)', ADVANCE = 'YES') &
' IN i: ', s1, '  S* (', Smin,  ', ', Smax, ')  |  nV = ', distanceInnerV
! READ(*,*) 

ELSE

WRITE(*, FMT = '(1X, A, I3, A, E9.3)', ADVANCE = 'YES') &
' IN i: ', s1, ' )  acceleration step |  nV = ', distanceInnerV
! READ(*,*) 

END IF

WRITE(*,*)
E1Stare = TE1Stare; AlphaStare = TAlphaStare; Sstare = TSstare; MStare = TMStare; ZiRawStare = TZiRawStare
V = TV
		
END DO				! CONTRACTION MAPPING STOPS HERE

CALL CPU_TIME(finish); finish = finish - begin

WRITE(*, '(1X, A, F12.4, A)') ' inner loop took ', finish, ' seconds '


END SUBROUTINE InvInnerE


!!	******************************************************************************************	!!
!!																								!!
!!																								!!
!!											InvOuter											!!
!!																								!!
!!																								!!
!!	******************************************************************************************	!!


SUBROUTINE InvOuter(alpha, delta, abeta, bbeta, masspoint, zibound, outputcost, thetam, thetan, &
					Atech, beta, Apref, storage, Nz, NS, NS0, Zvec, PI, knotsS, &
					rknotsV, knotsV, md, numelemV, cV, clmatV, LV, UV, dtauV, csV1, SVECV, SVECE1, &
					indicator, BetaKM, BetaSM, BetaP, BetaPQ, &
					precisionm, precisionS, precisionp, precisionq, simlength, zIsim, pbounds, &
					qbounds, KMbound, SMbound, Jmax, pvecnum, s1, PB, MUsim, Ssim, Ksim, &
					Qagg, Pagg, qmiss, Ydist, Mdist, NCdist, Cagg, Yagg, Xagg, Iagg, NKagg, &
					NMagg, Alphadist, AdjCost, distanceKinv, Qerror, Perror, &
					numberviolationsKM, numberviolationsSM, numberviolationsS, &
					numberviolationsQ, Sprodndist)

! Outer loop of the (S,s) inventory model.

INTEGER(ik):: Nz, NS, NS0, md, numelemV(4,md), simlength, s1, pvecnum, &
			  zIsim(simlength), Jmax, clmatV, t, zIval, jzf, samebins(Jmax+2), &
			  rknotsV(md), jS, jS0, Jloc(1), Kviolation, SMviolation, sviolation, &
			  qviolation, boundviolation, numberviolationsKM, numberviolationsSM, &
			  numberviolationsS,numberviolationsQ, si
			  	
REAL(rk):: alpha, delta, abeta, bbeta, masspoint, zibound, outputcost, thetam, thetan, &
		   gz, Atech, beta, Apref, storage, knotsS(NS0), knotsV(clmatV+2,md), &
		   SVECV(2,1), SVECE1(2,1), cV(4, numelemV(3,md), Nz), LV(1, NS0-3), UV(2, NS0-2), &
		   dtauV(NS0-1), precisionp, precisionq, PI(Nz, Nz), KMbound(2), SMbound(2), &
		   pbounds(2), qbounds(2), zval, gzval, KM, SM, Sval(1), Ksim(simlength+1), KMf, SMf, kval, &
		   NKagg(simlength), pickle, adjustors, Qagg(simlength), &
		   boop, boop2, Lp, Ytotal, Pagg(simlength), Xtotal, BetaKM(Nz,3), BetaSM(Nz,3), &
		   csV1(4,NS0-1), V1zf(Nz), V1(NS0), pl, pr, fl, fr, f, p, q, qmiss, precisionS, &
		   precisionm, csE1(4,NS0-1), Zvec(2,Nz), PB(simlength), BetaP(Nz,3), BetaPQ(Nz,3), &
		   Kf, Pvec(pvecnum), Fvec(pvecnum), Vf(1,1), Sf(1), Cagg(simlength), Xagg(simlength), &
		   Iagg(simlength), Yagg(simlength), Xsupply, AdjCost(simlength), NMtotal,  &
		   NMagg(simlength), distanceKinv, Qerror(simlength), Perror(simlength), ql, qr, &
		   qfl, qfr, Eal, Ear, Sstarl, Sstarr, E1error  !, &
		   !LE1S(1,1), LEa, Lq, LStest(1)
		   

		  	
REAL(rk)::	CYJ(Jmax+1), MYJ(Jmax+1), NYJ(Jmax+1)
REAL(rk)::	NCdist(simlength, Jmax+1), Ydist(simlength, Jmax+1), Mdist(simlength, Jmax+1), Sstar
REAL(rk)::	AlphaJ(Jmax+1)
REAL(rk)::	MUsim(simlength+1, Jmax+1), Ssim(simlength+1, Jmax+1), Alphadist(simlength, Jmax+1),  &
			ECostdist(simlength, Jmax+1), ECostJ(Jmax+1), z0val, &
			YTemp(Jmax+2), MTemp(Jmax+2), NTemp(Jmax+2), E1Temp(Jmax+2), E1(NS0), msol, nsol, cysol, E1S(1,1)

REAL(rk)::		SstarE1
REAL(rk)::		Sprodn(Jmax+1), Sfsim(Jmax+1), MUfsim(Jmax+1), Sprodndist(simlength, Jmax+1)


CHARACTER(20):: indicator
INTENT(IN):: alpha, delta, abeta, bbeta, masspoint, zibound, outputcost, thetam, thetan, &
			 Atech, beta, Apref, storage, Nz, NS, NS0, Zvec, PI, knotsS, &
			 rknotsV, knotsV, md, BetaP, BetaPQ, KMbound, numelemV, &
			 clmatV, LV, UV, dtauV, SVECV, SVECE1, indicator, BetaKM, BetaSM, precisionm, &
			 precisionS, cV, PB, zIsim, pvecnum, simlength, pbounds, qbounds, SMbound
INTENT(INOUT):: csV1, MUsim, Ssim, Ksim
INTENT(OUT):: Qagg, Pagg, qmiss, Ydist, Mdist, NCdist, Cagg, Yagg, Xagg, Iagg, NKagg, &
			  Alphadist, AdjCost
			  
INTENT(OUT)::	Sprodndist


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!							!!	
!!							!!
!!	Loop over Time Periods	!!
!!							!!
!!							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

numberviolationsKM = 0_ik
numberviolationsSM = 0_ik
numberviolationsS = 0_ik
numberviolationsQ = 0_ik

DO t = 1, simlength

! set the current state and forecasts of the future state.
zIval = ZIsim(t)
zval = Zvec(1,zIval)
gzval = Zvec(2,zIval)
KM = Ksim(t)
SM = DOT_PRODUCT(MUsim(t,:), Ssim(t,:))	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!														!
!	Expectations and Value Functions: The first step	!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

KMf = DOT_PRODUCT(BetaKM(zIval,1:3), (/1.0_rk, DLOG(KM), DLOG(SM)/)); KMf = DEXP(KMf)
SMf = DOT_PRODUCT(BetaSM(zIval,1:3), (/1.0_rk, DLOG(KM), DLOG(SM)/)); SMf = DEXP(SMf)


!! ********** IMPLEMENT GRID TRUNCATION ON FORECAST IF NECESSARY ************* !!
! This should never be necessary after model bounds have been found!
IF (KMf.GE.KMbound(2)) THEN 
 	WRITE(*,'(1X, A, F8.4, A, F8.4)', ADVANCE = 'YES') ' TRUNCATING KM FORECAST, (', KMf, '), as it VIOLATES UPPER Bound of ', KMbound(2)
 	WRITE(*,*)
	KMf = KMbound(2)
!
ELSEIF (KMf.LT.KMbound(1)) THEN
 	WRITE(*,'(1X, A, F8.4, A, F8.4)', ADVANCE = 'YES') ' TRUNCATING KM FORECAST, (', KMf, '), as it VIOLATES LOWER Bound of ', KMbound(1)
 	WRITE(*,*)
	KMf = KMbound(1)
END IF

IF (SMf.GE.SMbound(2)) THEN  
 	WRITE(*,'(1X, A, F8.4, A, F8.4)', ADVANCE = 'YES') ' TRUNCATING SM FORECAST, (', SMf, '), as it VIOLATES UPPER Bound of ', SMbound(2)
 	WRITE(*,*)
	SMf = SMbound(2)
!
ELSEIF (SMf.LT.SMbound(1)) THEN
 	WRITE(*,'(1X, A, F8.4, A, F8.4)', ADVANCE = 'YES') ' TRUNCATING SM FORECAST, (', SMf, '), as it VIOLATES LOWER Bound of ', SMbound(1)
 	WRITE(*,*)
	SMf = SMbound(1)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!																						    !!!	
!!!         WE FIX THE AGGREGATE STATE AND DERIVE V, THE UNIVARIATE FUNCTION ON S           !!!
!!!																						    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! convert multivariate splines to univariate splines, fit a univariate spline to 
! V[zi](S;KM,SM) having fixed (KM,SM) and zi.  This is the expectation of 
! V(S;z(j), KMf, SMf) given z(i), KM and SM.  It is needed in the subroutine 
! QfromP to first determine E1, and generate a univariate spline csE1.  
DO jS = 1,NS0
	Sval = knotsS(jS)
		
	DO jzf = 1,Nz
		V1zf(jzf) = FastSplinEval(md, numelemV(3,md), cV(:,:,jzf), rknotsV-1_ik, clmatV+2_ik, &
								  knotsV, (/Sval, KMf, SMf/), (/0_ik,0_ik, 0_ik/))
		V1zf(jzf) = PI(zIval, jzf)*V1zf(jzf)
	END DO
	
	V1(jS) = SUM(V1zf)
END DO

! Fit a spline to update polynomial approximations to E1 and V
CALL SPpp(V1, NS0 -2_ik, 1, LV, UV, dtauV, indicator, SVECV, csV1(:,:)) 

  	
!!!!!!!!!!!!!!!!!!!!!
!					!
!	Bisection on p	!
!					!
!!!!!!!!!!!!!!!!!!!!!

! Given the aggregate state (z(i), KM, SM) and thus the forecast of the future state, 
! (zf(j), KMf, SMf), we have the value of any inventory for any firm.  For any p, the 
! marginal utility of household consumption, and q, the relative price of the inter-
! mediate good, this will provide optimal choices including adjustment behavior.  
! We bisect for p, and within this for q, such that (1) the supply and demand of 
! intermediate goods is satisfied and (2) that p = MU(C).  

kval = KM; gz = Atech*gzval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!																						!
!	Our first attempt to bisect on p, locate two p' = f(p) such that f(pl)*f(pr)<0		!
!																						!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pl = pbounds(1)*PB(t)   !These are pbounds set up in the input file.  Here, trying
						!a state-dependent alternative that applies to the Baseline
						!pl = Pbound(zval, kval, alpha, Apref, delta, klow, precisionp, 1_ik)

pr = pbounds(2)*PB(t)	!These are pbounds set up in the input file.  Here, trying
						!a state-dependent alternative that applies to the Baseline
						!pr = Pbound(zval, kval, alpha, Apref, delta, khigh, precisionp, 0_ik)
Pvec = 0.0_rk; Fvec = 0.0_rk

CALL LINSPACE(pl, pr, pvecnum, Pvec)
pickle = 0.0_rk

DO jS0 = 1, pvecnum

	p = Pvec(jS0)

! Here we learn all we need to know about non-adjustors
! before calling up QfromP which will use, but not re-compute
! this information, as it does not depend on q.
YTemp = 0.0_rk; MTemp = 0.0_rk; NTemp = 0.0_rk; E1Temp = 0.0_rk; samebins = 0_ik
DO si = 1, Jmax+1

	IF (Ssim(t,si).GT.precisions) THEN
	
		CALL InvM(NS0, thetam, thetan, Atech, p, beta, precisionm, Apref, storage, 1.0_rk, knotsS, &
					 csV1(:,:), (/0.0_rk, Ssim(t,si)/), gzval, Ssim(t,si), msol, nsol, cysol, E1S(1,1))

	ELSE

		msol = Ssim(t,si); Sf(1) = 0.0_rk; CALL SPeval(csV1, knotsS, NS0 - 2_ik, 1, Sf, 1, 0, Vf)
		nsol = thetan*p*gz*(msol**thetam)/Apref
		nsol = nsol**(1.0_rk/(1.0_rk - thetan))
		cysol = gz*(msol**thetam)*(nsol**thetan)
		E1S(1,1) = p*cysol - Apref*nsol + beta*Vf(1,1)
		samebins(si+1) = 1_ik			! Svec next period for this bin will be zero

	END IF
																				
	YTemp(si+1) = cysol								! This vector will be the output of each nonadjusting
	MTemp(si+1) = msol								! group.  First entry still zero, because the Sstar
	NTemp(si+1) = nsol								! has not yet been determined, and this first entry
	E1Temp(si+1) = E1S(1,1)							! corresponds to adjustors.  Same holds for Mtemp, etc.	

END DO

! By the same token, the knot-point midperiod values [E1(knotsS)], and hence its spline
! should not depend on q, and hence the spline is computed just one time for any
! given p.

! This is just like the Inner Loop.  We use the univariate spline on V to create a grid of 
! E1 values, and then we fit a univariate spline to the E1, csE1.  Passed through QfromP, this
! spline is sent to EconomyQ, which must determine S*.   csV, a univariate spline on V, 
! must already be available at the start of this subroutine.

csE1 = 0.0_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!																!
!																!	
!	PRELIMINARY LOOK FOR BISECTION BOUNDS OF p - LP EQUATION	!
!																!
!																!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! retrieve E1 given p and fit a univariate spline to it
DO jS = 1, NS0
	Sval(1) = knotsS(jS)

	IF (Sval(1).GT.precisions) THEN		! AK NEW LINES 03.03.2003
	
		CALL InvM(NS0, thetam, thetan, Atech, p, beta, precisionm, Apref, storage, 1.0_rk, &
					 knotsS, csV1, (/0.0_rk, Sval(1)/), gzval, Sval(1),  msol, nsol, cysol, E1(jS))

	ELSE

		msol = Sval(1); Sf(1) = 0.0_rk; CALL SPeval(csV1, knotsS, NS0 - 2_ik, 1, Sf, 1, 0, Vf)
		nsol = thetan*p*gz*(msol**thetam)/Apref
		nsol = nsol**(1.0_rk/(1.0_rk - thetan))
		cysol = gz*(msol**thetam)*(nsol**thetan)
		E1(jS) = p*cysol - Apref*nsol + beta*Vf(1,1)
		
	END IF				
			
END DO

CALL SPpp(E1, NS0 -2_ik, 1, LV, UV, dtauV, indicator, SVECE1, csE1(:,:))

CALL QfromP(alpha, delta, abeta, bbeta, masspoint, zibound, outputcost, thetam, thetan, &
            Atech, beta, Apref, storage, NS0, Nz, Zvec, PI, zIval, &
			BetaP, BetaPQ, BetaSM, SM, knotsS, precisionq, csV1, &
			qbounds, KMbound, SMbound, Jmax, zval, gzval, kval, MUsim(t,:), Ssim(t,:), &
			samebins, p, precisionS, precisionm, csE1, Lp, Ytotal, NMtotal, &
			q, Xtotal, Kf, qmiss, adjustors, Sstar, ql, qr, qfl, qfr, Eal, Ear, Sstarl, &
			Sstarr, Sfsim, MUfsim, CYJ, MYJ, NYJ, AlphaJ, ECostJ, Sprodn, MTemp, &
			Ytemp, NTemp, E1Temp, SstarE1, E1error)

Fvec(jS0) = p - Lp

IF (jS0.GT.1_ik.AND.pickle.EQ.0.0_rk.AND.Fvec(jS0)*Fvec(MAX(jS0-1,1)).LT.0.0_rk) THEN
	pl = Pvec(jS0-1)
	pr = Pvec(jS0)
	fl = Fvec(jS0-1)
	fr = Fvec(jS0)
	pickle = 1.0_rk
	EXIT
ELSEIF(jS0.EQ.pvecnum.AND.pickle.EQ.0.0_rk) THEN
	pl = Pvec(1); pr = Pvec(pvecnum); fl = Fvec(1); fr = Fvec(pvecnum)
END IF

END DO

pickle = MAXVAL((Fvec(1:pvecnum - 1) - Fvec(2:pvecnum)))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!																						!	
!	If we managed to locate a negative and a positive value for f, then we can bisect	!
!				REAL BISECTION FOR p - Lp using the bounds found above					!
!																						!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


IF (fl*fr.GE.0.0_rk) THEN

 WRITE(*,'(1X, A, F8.4, A, F8.4, A, E12.4, A, E12.4)', ADVANCE = 'YES') ' pleft: ', pl, ' pright: ', pr, ' fleft: ',fl, ' fright: ', fr

ELSE 

	DO 	
		IF (pr - pl.LT.precisionp) THEN	
			EXIT
		END IF

		p = 0.5_rk*pl + 0.5_rk*pr

! Here we learn all we need to know about non-adjustors
! before calling up QfromP which will use, but not re-compute
! this information, as it does not depend on q.

YTemp = 0.0_rk; MTemp = 0.0_rk; NTemp = 0.0_rk; E1Temp = 0.0_rk; samebins = 0_ik

DO si = 1, Jmax+1

	IF (Ssim(t,si).GT.precisions) THEN	
	
		CALL InvM(NS0, thetam, thetan, Atech, p, beta, precisionm, Apref, storage, 1.0_rk, knotsS, &
					 csV1, (/0.0_rk, Ssim(t,si)/), gzval, Ssim(t,si), msol, nsol, cysol, E1S(1,1))

	ELSE
		
		msol = Ssim(t,si); Sf(1) = 0.0_rk; CALL SPeval(csV1, knotsS, NS0 - 2_ik, 1, Sf, 1, 0, Vf)
		nsol = thetan*p*gz*(msol**thetam)/Apref
		nsol = nsol**(1.0_rk/(1.0_rk - thetan))
		cysol = gz*(msol**thetam)*(nsol**thetan)
		E1S(1,1) = p*cysol - Apref*nsol + beta*Vf(1,1)
		samebins(si+1) = 1_ik		! This Svec tomorrow for this bin will be zero
	
	END IF
																				
	YTemp(si+1) = cysol								! This vector will be the output of each nonadjusting
	MTemp(si+1) = msol								! group.  First entry still zero, because the Sstar
	NTemp(si+1) = nsol								! has not yet been determined, and this first entry
	E1Temp(si+1) = E1S(1,1)							! corresponds to adjustors.  Same holds for Mtemp, etc.	

END DO

csE1 = 0.0_rk

! retrieve E1 given p and fit a univariate spline to it
DO jS = 1, NS0
	Sval(1) = knotsS(jS)

	IF (Sval(1).GT.precisions) THEN		
		
		CALL InvM(NS0, thetam, thetan, Atech, p, beta, precisionm, Apref, storage, 1.0_rk, &
				  knotsS, csV1, (/0.0_rk, Sval(1)/), gzval, Sval(1),  msol, nsol, cysol, E1(jS))
	ELSE

		msol = Sval(1); Sf(1) = 0.0_rk; CALL SPeval(csV1, knotsS, NS0 - 2_ik, 1, Sf, 1, 0, Vf)
		nsol = thetan*p*gz*(msol**thetam)/Apref
		nsol = nsol**(1.0_rk/(1.0_rk - thetan))
		cysol = gz*(msol**thetam)*(nsol**thetan)
		E1(jS) = p*cysol - Apref*nsol + beta*Vf(1,1)
		
	END IF				
			
END DO

CALL SPpp(E1, NS0 -2_ik, 1, LV, UV, dtauV, indicator, SVECE1, csE1(:,:))

CALL QfromP(alpha, delta, abeta, bbeta, masspoint, zibound, outputcost, thetam, thetan, &
		    Atech, beta, Apref, storage, NS0, Nz, Zvec, PI, zIval, &
			BetaP, BetaPQ, BetaSM, SM, knotsS, precisionq, csV1, &
			qbounds, KMbound, SMbound, Jmax, zval, gzval, kval, MUsim(t,:), Ssim(t,:), &
			samebins, p, precisionS, precisionm, csE1, Lp, Ytotal, NMtotal, &
			q, Xtotal, Kf, qmiss, adjustors, Sstar, ql, qr, qfl, qfr, Eal, Ear, Sstarl, &
			Sstarr, Sfsim, MUfsim, CYJ, MYJ, NYJ, AlphaJ, ECostJ, Sprodn, MTemp, &
			Ytemp, NTemp, E1Temp, SstarE1, E1error)

f = p - Lp

IF (f*fl.GT.0.0_rk) THEN
	pl = p; fl = f
ELSE
	pr = p; fr = f
END IF

END DO

END IF
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							!
!	Recover the economy		!
!							!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We now have equilibrium p and q 
Ksim(t+1) = Kf

! Additions to Stored Simulation Results:
Sprodndist(t,:) = Sprodn; Ssim(t+1,:) = Sfsim; MUsim(t+1,:) = MUfsim
Ydist(t,:) = CYJ; Mdist(t,:) = MYJ; NCdist(t,:) = NYJ
Qagg(t) = q
Pagg(t) = p

! Moved definition of Xsupply up in code, specifically so as to set Xagg(t) as Xsupply (rather than demand)
z0val = zval

Xsupply = p*(1.0_rk - alpha)/Apref
Xsupply = Xsupply**((1.0_rk - alpha)/alpha)
Xsupply = (z0val**(1.0_rk/alpha))*Xsupply*kval
Xsupply = Xsupply*(q**((1.0_rk - alpha)/alpha))

Xagg(t) = Xsupply  !rather than Xtotal

Iagg(t) = Ksim(t+1) - (1.0_rk - delta)*Ksim(t)
Alphadist(t,:) = AlphaJ
ECostdist(t,:) = ECostJ

NKagg(t) = p*q*z0val*(1.0_rk - alpha)/Apref; NKagg(t) = (NKagg(t)**(1.0_rk/alpha))*kval
NMagg(t) = NMtotal

AdjCost(t) = DOT_PRODUCT(ECostdist(t,:), MUsim(t,:))
Yagg(t) = Ytotal
Cagg(t) = Ytotal - Ksim(t+1) + (1.0_rk - delta)*Ksim(t)

Jloc = MAXLOC(AlphaJ)							! Must be a single dimensional array
boop = DOT_PRODUCT(MUsim(t+1,:), Ssim(t+1,:))
boop2 = DOT_PRODUCT(MUsim(t+1,1:Jmax+1), MYJ) 

		
WRITE(*,*) ' ' 
WRITE(*, '(1X, I3, A, I5, A, I2, A, F5.2, 2(A, F6.3), A, 2(F7.3), A, 2(F7.3), A, I2)', &
ADVANCE = 'YES' ) s1, ' | ', t, ' z', zIval, ' p ', Pagg(t), ' /p0 ', Pagg(t)/PB(t), &
' q ', q, ' (kf, SMf) ', Ksim(t+1), boop, '  S*i ', Sstar, SstarE1, '  J ', Jloc(1)

WRITE(*,'(1X, 2(A,F6.4), (A,2F8.4), (A, 2E9.2), 2(A, E9.2))', ADVANCE = 'YES') '  XS: ', &
		Xsupply, ' XD: ', Xtotal, ' (Eal, Ear): ', Eal, Ear, ' (qfl, qfr): ', qfl, qfr, '  qE  ', qmiss, ' dk ', distanceKinv

WRITE(*,*) ' '
WRITE(*,'(1X, 4(A, F10.6), A)') ' E1error = ', E1error, '  Adjust = ',adjustors, '  Y: ', &
        Ytotal, ' M: ', boop2, ' Alpha(t), S1(t), S(t+1) and MU(t+1):'

WRITE(*,'(1X,A)', ADVANCE = 'NO') '                 '

DO jS = 1, Jmax-1
	WRITE(*,'(1X,F8.3)', ADVANCE = 'NO') ALPHAJ(jS)	! adjustment fractions
END DO


WRITE(*,'(1X,A)', ADVANCE = 'YES') '  '
DO jS = 1, Jmax+1
	WRITE(*,'(1X,F8.3)', ADVANCE = 'NO') Sprodn(jS)		! production time S1(t)
END DO


WRITE(*,'(1X,A)', ADVANCE = 'YES') '  '
DO jS = 1, Jmax+1
	WRITE(*,'(1X,F8.3)', ADVANCE = 'NO') Ssim(t+1,jS)	! S(t+1) directly associated with S1(t) above
END DO


WRITE(*,'(1X,A)', ADVANCE = 'YES') '  '
DO jS = 1, Jmax+1
	WRITE(*,'(1X,F8.3)', ADVANCE = 'NO') MUsim(t+1,jS)	! MU(t+1) directly associated with S(t+1) above
END DO

WRITE(*,'(1X,A, F12.4)', ADVANCE = 'YES') '  ', SUM(MUsim(t+1,:))

Qerror(t) = qmiss
Perror(t) = f

!!!!!!!! REPORTING BOUND VIOLATIONS SECTION !!!!!!!!!!
boundviolation = 0_ik
Kviolation = 0_ik
SMviolation = 0_ik
sviolation = 0_ik
qviolation = 0_ik

! I put in a '5' for any upper knot violation, a '1' for any lower knot violation
IF (Ksim(t+1).GT.KMbound(2).OR.Ksim(t+1).LT.KMbound(1)) THEN
	boundviolation = boundviolation + 1
	IF (Ksim(t+1).GT.KMbound(2)) THEN; Kviolation = 5_ik; END IF
	IF (Ksim(t+1).LT.KMbound(1)) THEN; Kviolation = 1_ik; END IF
	numberviolationsKM = numberviolationsKM + 1
END IF

!SM index = 3
IF (boop.GT.SMbound(2).OR.boop.LT.SMbound(1)) THEN
	boundviolation = boundviolation + 1
	IF (boop.GT.SMbound(2)) THEN; SMviolation = 5_ik; END IF
	IF (boop.LT.SMbound(1)) THEN; SMviolation = 1_ik; END IF
	numberviolationsSM = numberviolationsSM + 1
END IF
	

IF (Sstar.GE.knotsS(NS0).OR.Sstar.LE.knotsS(1)) THEN
	boundviolation = boundviolation + 1
	IF (Sstar.GE.knotsS(NS0)) THEN; sviolation = 5_ik; END IF
	IF (Sstar.LE.knotsS(1)) THEN; sviolation =  1_ik; END IF	
	numberviolationsS = numberviolationsS + 1
END IF

IF (q.GE.qbounds(2)-precisionp.OR.q.LE.qbounds(1)+precisionp) THEN
	boundviolation = boundviolation + 1
	IF (q.GE.qbounds(2)-precisionp) THEN; qviolation = 5_ik; END IF
	IF (q.LE.qbounds(1)+precisionp) THEN; qviolation =  1_ik; END IF	
	numberviolationsQ = numberviolationsQ + 1
END IF

IF (boundviolation.GT.0_ik) THEN
	WRITE(*, '(1X,A, 4(I4), 1X,A)', ADVANCE = 'YES') &
	 '       Violations on K,SM,S*,q:', Kviolation, SMviolation, sviolation, qviolation, &
	 '              (violations: 5=upper, 1=lower)'		           

WRITE(*,*) ' '

END IF
!!!!!!!! BOUND VIOLATION SECTION !!!!!!!!!!


END DO			! THE NEXT TIME PERIOD


END SUBROUTINE InvOuter


!!	******************************************************************************************	!!
!!																								!!
!!																								!!
!!											QfromP												!!
!!																								!!
!!																								!!
!!	******************************************************************************************	!!

!!	QfromP: q = Q(p) such that XD = XS	

SUBROUTINE QfromP(alpha, delta, abeta, bbeta, masspoint, zibound, ouputcost, thetam, thetan, &
				  Atech, beta, Apref, storage, NS0, Nz, Zvec, PI, zi, &
				  BetaP, BetaPQ, BetaSM, SM, knotsS, precisionq, &
				  csV, qbounds, KMbound, SMbound, Jmax, zval, gzval, kval, muvec, Sivec, &
				  samebins, p, precisions, precisionm, csE1, Lp, Ytotal, &
				  NMtotal, q, Xtotal, Kf, qmiss, adjustors, Sstar, ql, qr, fl, fr, Eal, &
				  Ear, Sstarl, Sstarr, Sfsim, MUfsim, CYJ, MYJ, NYJ, AlphaJ, ECostJ, &
				  Sprodn, MTemp, Ytemp, NTemp, E1Temp, SstarE1, E1error)


INTEGER(ik):: NS0, Jmax, s1, Nz, zi, qvecnum, jq, jabba, samebins(Jmax+2), js0
REAL(rk):: alpha, delta, thetam, thetan, Atech, beta, Apref, storage, knotsS(NS0), &
		   precisionq, csV(4,NS0-1), qbounds(2), KMbound(2), SMbound(2), abeta, bbeta, &
		   masspoint, zibound, ouputcost, zval, gzval, csE1(4,NS0-1), p, precisions, precisionm, kval, Lp, &
		   Ytotal, SstarE1, NMtotal, q, Kf, BetaP(Nz,3), BetaPQ(Nz,3), BetaSM(Nz,3), Cp, &
		   Zvec(2,Nz), PI(Nz,Nz), SM, Xtotal, SHIFTJ(Jmax-1), qmiss, termqsearch, SMfp, MYJ(Jmax + 1), &
		   CYJ(Jmax + 1), NYJ(Jmax + 1), adjustors, pickle, ql, fl, qr, fr, f, Ea, Eal, Ear, &
		   Sstarl, Sstarr, myjl, myjr, cyjl, cyjr, nyjl, nyjr, Sfsim(Jmax + 1), MUfsim(Jmax + 1), &
		   MTemp(Jmax + 2), Ytemp(Jmax + 2), NTemp(Jmax + 2), E1Temp(Jmax + 2), E1error, zerobinmeasure, &
		   Z0vec(Nz)

REAL(rk), ALLOCATABLE:: Qvec(:), FQvec(:)

REAL(rk)::	Sivec(Jmax+1), MUvec(Jmax+1), AlphaJ(Jmax+1), ECostJ(Jmax+1), Sstar

REAL(rk)::		Sprodn(Jmax+1)

INTENT(IN):: alpha, delta, abeta, bbeta, masspoint, zibound, ouputcost, thetam, thetan, &
			 Atech, beta, Apref, NS0, knotsS, precisionq, csV, qbounds, Jmax, &
			 zval, gzval, kval, muvec, Sivec, p, precisions, Nz, Zvec, PI, zi, BetaP, BetaPQ, BetaSM, &
			 SM, csE1, MTemp, Ytemp, NTemp, E1Temp, storage, samebins

INTENT(OUT):: Ytotal, NMtotal, q, Xtotal, Kf, qmiss, Lp, adjustors, Sstar, ql, qr, fl, fr, Eal, Ear, &
			  Sstarl, Sstarr, Sfsim, MUfsim, MYJ, NYJ, CYJ, AlphaJ, ECostJ, Sprodn, SstarE1, E1error


s1 = 0_ik
Eal = 0.0_rk; Ear = 0.0_rk; Sstarl = 0.0_rk; Sstarr = 0.0_rk

qvecnum = 5_ik

! numerical term that is the constant part of the right hand side of the equation in step 2.  
! Z(1,:) = z values (intermediate goods baseline shock), Z(2,:) = gz values (final goods shock)

DO jabba = 1, Nz
    Z0vec(jabba) = Zvec(1,jabba)
END DO

termqsearch = p*(1.0_rk - alpha)/Apref
termqsearch = termqsearch**((1.0_rk - alpha)/alpha)
termqsearch = (zval**(1.0_rk/alpha))*termqsearch*kval

!! ************************************ !!
!!										!!
!!	Search for a bisection in q given p	!!
!!										!!
!! ************************************ !!

ALLOCATE(Qvec(qvecnum), FQvec(qvecnum))

Qvec = 0.0_rk; FQvec = 0.0_rk

CALL LINSPACE(qbounds(1), qbounds(2), qvecnum, Qvec)

jabba = 0; pickle = 0.0_rk

DO jq = 1, qvecnum	

	q = Qvec(jq)

	CALL DynEconomyQ(NS0, p, q, precisions, precisionm, knotsS, csE1, abeta, bbeta, masspoint, zibound, &
			         ouputcost, Apref, thetam, thetan, gzval, Atech, beta, storage, csV, &
				     Jmax, MUvec, Sivec, AlphaJ, Sstar, SstarE1, Ea, MYJ, CYJ, NYJ, adjustors, Ytotal, &
				     NMtotal, Xtotal, ECostJ, MTemp, Ytemp, NTemp, E1Temp, E1error)

	qmiss = Xtotal - termqsearch*(q**((1.0_rk - alpha)/alpha)); FQvec(jq) = qmiss

	IF (qmiss.LT.0.0_rk) THEN
		Eal = Ea
		Sstarl = SstarE1
	ELSE
		Ear = Ea
		Sstarr = SstarE1
	END IF

	IF (jq.GT.1_ik.AND.pickle.EQ.0.0_rk.AND.FQvec(jq)*FQvec(MAX(jq-1,1)).LT.0.0_rk) THEN
		ql = Qvec(jq-1)
		fl = FQvec(jq-1)
		qr = Qvec(jq)
		fr = FQvec(jq)
		pickle = 1.0_rk
		EXIT
	END IF
	
END DO

IF (pickle.EQ.0.0_rk) THEN
	
	WRITE(*,'(1X,A,F8.4,A)') ' Unable to bisect q given p = ', p, ' qbounds must be extended.  Printing FQvec ' 	 
	PRINT*, FQvec
	READ(*,*)

END IF

DEALLOCATE(Qvec, FQvec)


!!!!!!!!!!!!!!!!!!!!!!!!!
!						!		
!		BISECTION		!
!						!
!!!!!!!!!!!!!!!!!!!!!!!!!

DO 
	
	IF (qr - ql.LT.precisionq) THEN				
		EXIT
	END IF
	
	q = 0.5_rk*ql + 0.5_rk*qr

	CALL DynEconomyQ(NS0, p, q, precisions, precisionm, knotsS, csE1, abeta, bbeta, masspoint, zibound, &
			         ouputcost, Apref, thetam, thetan, gzval, Atech, beta, storage, csV, &
				     Jmax, MUvec, Sivec, AlphaJ, Sstar, SstarE1, Ea, MYJ, CYJ, NYJ, adjustors, Ytotal, &
				     NMtotal, Xtotal, ECostJ, MTemp, Ytemp, NTemp, E1Temp, E1error)
	
	! XS(q) = termqsearch*(q**(1 - alpha)/alpha) and thus f(q) = XD(q) - XS(q), the excess demand for X at q.
	qmiss = Xtotal - termqsearch*(q**((1.0_rk - alpha)/alpha))
	f = qmiss

	IF (f*fl.GT.0.0_rk) THEN
		ql = q; fl = f; Sstarl = Sstar; Eal = Ea; myjl = MYJ(1); cyjl = CYJ(1); nyjl = NYJ(1)
	ELSE
		qr = q; fr = f; Sstarr = Sstar; Ear = Ea; myjr = MYJ(1); cyjr = CYJ(1); nyjr = NYJ(1)
	END IF

END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!												!
!												!
!		RETRIEVE EQUILIBRIUM DISTRIBUTION		!
!												!
!												!			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Sprodn = 0.0_rk; Sfsim = 0.0_rk; MUfsim = 0.0_rk; SHIFTJ = 0.0_rk

! MYJ, CYJ and NYJ are unaltered from the call to EconomyQ
Sfsim(1) = Sstar - MYJ(1); MUfsim(1) = DOT_PRODUCT(AlphaJ,MUvec)
	
js0 = 2_ik; zerobinmeasure = 0.0_rk

DO jq = 2, Jmax + 1
			
	IF (samebins(jq).NE.1_ik) THEN		! This is not a bin that implies [m(j) = s(j) hence s'(j) = 0.0]

		MUfsim(js0) = (1.0_rk - AlphaJ(jq-1))*MUvec(jq-1)
		Sfsim(js0) = Sivec(jq-1) - MYJ(jq)
			
		js0 = js0 + 1_ik

	ELSE							! This is a bin that has s'(j) = 0.0

		zerobinmeasure = zerobinmeasure + (1.0_rk - AlphaJ(jq-1))*MUvec(jq-1)
		
	END IF 

END DO

IF (js0.LT.Jmax + 1_ik) THEN		! aggregate all bins with s'(j) = 0.0
	MUfsim(js0) = zerobinmeasure
END IF

Sprodn = (/Sstar, Sivec(1:Jmax)/)

! Note that Ytotal, from EconomyQ in InvToolShed, is already net of any relevant output-denominated 
! adjustment costs.  Hence given total production based on (p,q) we have net output of final 
! goods, Ytotal.  

Cp = 1.0_rk/p								  ! C that's implied by P
Kf = (1.0_rk - delta)*kval + (Ytotal - Cp)    !True Kf (rather than forecast) for next period

! If we put in 0.0_rk, ValueKf uses the forecasted SMf, rather than actual.
SMfp = DOT_PRODUCT(MUfsim, Sfsim)
Lp = ValueKf(Nz, Z0vec, PI, zi, alpha, delta, Apref, BetaP, beta, Kf, BetaPQ, BetaSM, kval, SM, KMbound, SMbound, SMfp)


END SUBROUTINE QfromP


!!	******************************************************************************************	!!
!!																								!!
!!																								!!
!!											ValueKf						 					    !!
!!																								!!
!!																								!!
!!	******************************************************************************************	!!

!_____________________________
!    
!     ValueKf PURE FUNCTION
!_____________________________

PURE FUNCTION ValueKf(Nz, Z0vec, PI, zi, alpha, delta, Apref, BetaP, beta, Kf, BetaPQ, BetaSM, kval, SM, KMbound, SMbound, SMf)

INTEGER(ik):: Nz, zi, zfi
REAL(rk):: Z0vec(Nz), PI(Nz, Nz), alpha, delta, Apref, BetaP(Nz,3), BetaPQ(Nz,3), BetaSM(Nz,3), &
		   beta, Kf, Pf(Nz), term(Nz), term0, term1, ValueKf, SM, SMf, PQf(Nz), kval, &
		   KMbound(2), SMbound(2), Kf0, SMf0, SMf1
INTENT(IN):: Nz, Z0vec, PI, zi, alpha, delta, Apref, BetaP, beta, Kf, BetaPQ, BetaSM, &
             kval, SM, KMbound, SMbound, SMf

term0 = (1.0_rk - alpha)/Apref  
term1 = (1.0_rk - alpha)/alpha
term0 = term0**term1
term0 = alpha*term0

IF (SMf.EQ.0.0_rk) THEN		! If we DID NOT put in the real implied SMf, but rather want it forecast-determined)
	SMf1 = DOT_PRODUCT(BetaSM(zi,1:3), (/1.0_rk, DLOG(kval), DLOG(SM)/)); SMf1 = DEXP(SMf)
ELSE
	SMf1 = SMf
END IF

	!! ********** IMPLEMENT GRID TRUNCATION IF NECESSARY ************* !!
	IF (Kf.GE.KMbound(2)) THEN 
		Kf0 = KMbound(2)
	ELSEIF (Kf.LT.KMbound(1)) THEN
		Kf0 = KMbound(1)
	ELSE 
		Kf0 = Kf
	END IF

	IF (SMf1.GE.SMbound(2)) THEN  
		SMf0 = SMbound(2)		

	ELSEIF (SMf1.LT.SMbound(1)) THEN
		 SMf0 = SMbound(1)
	ELSE
		SMf0 = SMf1
	END IF
	!! ********** IMPLEMENT GRID TRUNCATION IN NECESSARY ************* !!

DO zfi = 1, Nz

    Pf(zfi) = DOT_PRODUCT(BetaP(zfi,:), (/1.0_rk, DLOG(Kf0), DLOG(SMf0)/)); Pf(zfi) = DEXP(Pf(zfi))
	PQf(zfi) = DOT_PRODUCT(BetaPQ(zfi,:), (/1.0_rk, DLOG(Kf0), DLOG(SMf0)/)); PQf(zfi) = DEXP(PQf(zfi))
	
	term(zfi) = Pf(zfi)*(1.0_rk - delta) + term0*((PQf(zfi)*Z0vec(zfi))**(1.0_rk/alpha))
	term(zfi) = beta*PI(zi,zfi)*term(zfi)
END DO

ValueKf = SUM(term)


END FUNCTION ValueKf


!!	******************************************************************************************	!!
!!																								!!
!!																								!!
!!			Economy Q Subroutine: Determine Final Goods Firms' behavior given q				!!
!!																								!!
!!																								!!
!!	******************************************************************************************	!!


SUBROUTINE DynEconomyQ(NS, p, q, precisions, precisionm, knotsS, csE1, abeta, bbeta, masspoint, zibound, &
					outputcost, Apref, thetam, thetan, gzval, Atech, beta, storage, &
					csV, Jmax, MUvec, Svec, AlphaJ, Sstar, SstarE1, Ea, MYJ, CYJ, NYJ, adjustors, &
					Ytotal, NMtotal, Xtotal, ECostJ, MTemp, Ytemp, NTemp, E1Temp, E1error)

! This subroutine does more or less everything in the consumption-goods
! sector.  It takes as given p and q. It requires univariate function 
! csV for V that already reflects expectations of (zj, KMf, SMf) across 
! firms.  Further it requires a univariate function that captures 
! E1, csE1, the middle of period value, before production, of having some 
! mid-period level of inventory holdings, S1.  

INTEGER(ik):: NS, Jmax, si, jS
 
REAL(rk)::	p, q, precisions, precisionm, knotsS(NS), csE1(4,NS-1), Apref, thetam, &
			thetan, gzval, Atech, beta, csV(4,NS-1), adjustors, Ytotal, &
			abeta, bbeta, masspoint, zibound, outputcost, Ea, cysol, nsol, msol, &
			E1S(1,1), zitilde, Xtotal, costterm, NMtotal, E1error, storage, Ystar, &
			Mstar, Nstar		

REAL(rk)::	Sstar, MTemp(Jmax+2), Ytemp(Jmax+2), NTemp(Jmax+2), E1Temp(Jmax+2), &
			MUvec(Jmax+1), Svec(Jmax+1), AlphaJ(Jmax+1), ZIJraw(Jmax+1), ZIJ(Jmax+1), ECostJ(Jmax+1), &
			MYJ(Jmax+1), CYJ(Jmax+1), NYJ(Jmax+1), Sval(1)

INTEGER(ik):: Lnumlook, PlaceManualSstarE1(1)
REAL(rk)::	Lfraction, SstarLowE1, SstarUpE1, SstarE1
REAL(rk), ALLOCATABLE:: LvecS(:), LObjE1(:), LobjV(:), L2E1(:,:)
INTENT(IN):: NS, p, q, precisions, knotsS, csE1, Apref, thetam, thetan, &
			 Atech, beta, gzval, csV, Jmax, MUvec, Svec, abeta, bbeta, masspoint, zibound, &
			 outputcost, MTemp, Ytemp, NTemp, E1Temp

INTENT(OUT)::	AlphaJ, Sstar, Ea, MYJ, CYJ, NYJ, adjustors, Ytotal, NMtotal, Xtotal, ECostJ, SstarE1, E1error


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               																		!
!	What Consumption Goods Producers (Inventory-Holding Firms) Do at this P and Q?		!
!																		                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialize vectors
AlphaJ = 0.0_rk; MYJ = 0.0_rk; CYJ = 0.0_rk; NYJ = 0.0_rk; ECostJ = 0.0_rk

! Given our fit to E1, which comes from elsewhere as it is insensitive to q, we 
! determine S* and the cysol and msol of thoese who adjust to this level of inventory.
Lnumlook = 500_ik
ALLOCATE(LvecS(0:Lnumlook), LObjE1(0:Lnumlook), LObjV(0:Lnumlook), L2E1(2,0:Lnumlook))

DO jS = 0, Lnumlook
	Lfraction = DBLE(Lnumlook - jS)/DBLE(Lnumlook)
	LvecS(jS) = knotsS(1)*Lfraction + (1.0_rk - Lfraction)*knotsS(NS)
END DO
		
CALL SPeval(csE1, knotsS, NS-2_ik, 1, LvecS, Lnumlook+1, 0, LObjE1) 

LObjE1 = -(p*q)*LvecS + LObjE1; LObjE1 = 10000.00_rk*LObjE1
L2E1(1,:) = LvecS; L2E1(2,:) = LObjE1

!  Pick out SstarLowBd and SstarUpBd to send into the InvSstarSUP search.
PlaceManualSstarE1 = MAXLOC(LObjE1)

SstarE1 = LvecS(PlaceManualSstarE1(1)-1)

IF (PlaceManualSstarE1(1).GE.2_ik) THEN
	SstarLowE1 = LvecS(PlaceManualSstarE1(1) - 2_ik)
ELSE
	SstarLowE1 = LvecS(0)
END IF

IF (PlaceManualSstarE1(1).LE.Lnumlook - 1_ik) THEN
	SstarUpE1 = LvecS(PlaceManualSstarE1(1))
ELSE
	SstarUpE1 = LvecS(Lnumlook)
END IF


DEALLOCATE(LvecS, LObjE1, L2E1, LObjV)

! Call InvSstarSUP with the two last arguments being the lower and upper bound for the region
! of S1 that should be searched -- as isolated above.
CALL InvSstar(NS, p*q, precisions, knotsS, csE1, (/SstarLowE1, SstarUpE1/), Sval(1), Ea)
Sstar = Sval(1)

CALL InvM(NS, thetam, thetan, Atech, p, beta, precisionm, Apref, storage, 1.0_rk, knotsS, &
          csV, (/knotsS(1), Sstar/), gzval, Sstar, msol, nsol, cysol, E1S(1,1))

Ystar = cysol; Mstar = msol; Nstar = nsol

! By definition Ea = -p*q*Sstar + E1(Sstar)
IF (DABS(-p*q*Sstar + E1S(1,1) - Ea).GE.0.01_rk) THEN
	WRITE(*,*)
	WRITE(*,'(1X, A, F8.4)', ADVANCE = 'YES') ' Outer: Line #1190: Error -p*q*Sstar+E1(Sstar)-Ea: ', -p*q*Sstar + E1S - Ea
END IF

E1error = Ea - (-p*q*Sstar + E1S(1,1))

costterm = (1.0_rk - outputcost)*Apref + outputcost*p


!***********************************************************************************************!
! A NOTE ON THE PASSAGE BELOW, AND LAYOUT OF ALPHA, SVEC AND MUVEC RELATIVE TO CYJ, MYJ, NYJ	!
!																								!
!	At the start of the period, there are 1,...,Jmax+1 groups, hence the size of Alpha, and the	!
!	date t start-of-pd sizes of Svec and MUvec.													!
!																								!
!	There are Jmax+2 potential Svalues for production time: first, Sstar, and	                !
!	then the Svec values (which number Jmax+1).  Hence sizes of temporary prodn-time vectors	!
!	Ytemp, MTemp, and NTemp are all Jmax+2.														!
!																								!
!	However, Alpha(Jmax) = Alpha(Jmax+1) = 1.	                								!
!	This means that the only *relevant* elements of	[Sstar Svec], and the (size		            !
!	Jmax+2) vectors Ytemp, Mtemp, NTemp, SfutTemp will be Jmax+1 elments.  So, once we			!
!	have obtained these vectors, we will end up truncating them to send out only (Jmax+1)		!
!	relevant entries from each.  This is how the CYJ, MYJ, NYJ, Sprodn, Sfuture, MUfuture that	!
!	are sent out of this routine all have size(Jmax+1).											!
!																								!
!***********************************************************************************************!

DO si = 1, Jmax+1		! Find out what fraction of each group adjusts. E1S is independent of q.

	E1S(1,1) = E1Temp(si+1)
	zitilde = (p*q*Svec(si) - E1S(1,1) + Ea)/costterm	! This is equation (23) on page 6.
	ZIJraw(si) = zitilde							

	IF (si.LT.Jmax) THEN
		ZIJ(si) = DMIN1(DMAX1(0.0_rk, zitilde),zibound)
	ELSE
		ZIJ(si) = zibound							! alpha is 1 for Jmax and Jmax+1 groups
	END IF											! This must not bind.

END DO

! Use of 5th argument in Cost Function returns probability
DO si = 1, Jmax+1
	AlphaJ(si) = CostFunction(ZIJ(si), abeta, bbeta, masspoint, zibound, 0_ik)
	ECostJ(si) = CostFunction(ZIJ(si), abeta, bbeta, masspoint, zibound)
END DO													! filled in once 'h' and 'adjustors' are known below.

Xtotal = SUM((Sstar - Svec)*AlphaJ*MUvec)				! Demand from High Sstar adjustors

! Adjusting firms reach S* before production and produce CYJ(1). Those with S(j) that do not adjust 
! produce CYJ(j+1).

adjustors = SUM(AlphaJ*MUvec)

Ytotal = SUM((/adjustors, (1.0_rk - AlphaJ(1:Jmax+1))*MUvec(1:Jmax+1)/)*(/Ystar, Ytemp(2:Jmax+2)/)) 
Ytotal = Ytotal - outputcost*DOT_PRODUCT(ECostJ, MUvec)

NMtotal = SUM((/adjustors, (1.0_rk - AlphaJ(1:Jmax+1))*MUvec(1:Jmax+1)/)*(/Nstar, NTemp(2:Jmax+2)/))

MYJ = (/Mstar, MTemp(2:Jmax+1)/)	
CYJ = (/Ystar, YTemp(2:Jmax+1)/)
NYJ = (/Nstar, NTemp(2:Jmax+1)/)

END SUBROUTINE DynEconomyQ

END MODULE InvWorkShededit