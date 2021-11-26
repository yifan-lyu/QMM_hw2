 PROGRAM SteadyState

! Steady state for the inventory model of Khan and Thomas (2007).

! This program computes the steady state of a model with a nontrivial distribution of final 
! goods firms that hold different stocks of inventories.  The steady state is solved using 
! univariate piecewise polynomial splines.  

USE KindSet
use ppsplinefit3edit

IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                               !
!   Interface Block for external subroutines    !
!                                               !    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

interface

subroutine InvSstar(rknotsS, pq, precisionc, knotsS, csE1, Sbound, starscalar, Ea)

USE KINDSET
INTEGER(ik):: rknotsS
REAL(rk):: pq, precisionc, knotsS(rknotsS), Sbound(2), csE1(4,rknotsS-1), &
		   Ea, starscalar 

end subroutine invsstar

subroutine InvM(rknotsS, thetam, thetan, Atech, p, beta, precisiond, Apref, storage, storageh, knotsS, &
                csS, Sbound, zval, Sval, msol, nsol, cysol, E1sol)
USE KINDSET

INTEGER(ik):: rknotsS
REAL(rk):: thetam, thetan, Atech, p, beta, precisiond, Apref, storage, storageh, zval, &
		   Sval, knotsS(rknotsS), Sbound(2), csS(4,rknotsS-1), msol, nsol, &
		   E1sol, cysol 

end subroutine invm

end interface


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                   !
!   Type Declarations for Program   !
!                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER(ik):: NS, rknotsS, rknotsSpoly, rknotsSint, Jmax, maxJ, si, s1
INTEGER(ik), ALLOCATABLE:: samebins(:)
REAL(rk):: alpha, thetam, thetan, Atech, delta, Apref, Csupply, &
		   beta, storage, storage1, precisionv, precisiond, precisionc, SVECV(2,1), &
		   SVECE1(2,1), Sbound(2), kstar, p, q, Nkstar, M, Xnumadjust, abeta, bbeta, &
		   zibound, masspoint, zbar, Estar, Sstar, plow, phigh, Tplow, Tphigh, Tp0, distance, &
		   begintime, endingtime, outputcost, GDP, realrate, cyratio, kyratio, myratio, adjustors, &
		   Totalhours, plsuggest, phsuggest, logbound, Nm, AdCost
		   
REAL(rk), ALLOCATABLE:: knotsS(:), V(:), E1(:), csV(:,:), csE1(:,:), LW(:,:), UW(:,:), ECostJ(:), RWJ(:), &
						dtauW(:), SJ(:), AlphaJ(:), ThetaJ(:), MYJ(:), E1J(:), ZIJraw(:), NYJ(:), CYJ(:)
CHARACTER(30):: indicator, datafile, resultfile, startfile, storeq, storeh

!!!!!!!!!!!!!!!!!!!!!!!!!
!						!
!	Obtain Parameters	!
!						!
!!!!!!!!!!!!!!!!!!!!!!!!!

CALL SETUP(alpha, thetam, thetan, Atech, delta, abeta, bbeta, zibound, masspoint, &
		   outputcost, zbar, Apref, beta, storage, precisionv, precisiond, precisionc, NS, Sbound, &
		   Jmax, datafile, resultfile, startfile, plsuggest, phsuggest, logbound, indicator, SVECV, SVECE1) 

! default is original formulation of output-weighted storage cost, not as a value of inventory stock
storeq = 'y'; 
storeh = 'n';

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								!
!	 Univariate Spline Setup	!
!								!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rknotsS = NS
rknotsSpoly = rknotsS - 1_ik
rknotsSint = rknotsS - 2_ik

ALLOCATE(knotsS(rknotsS), LW(1, rknotsSint-1), UW(2, rknotsSint), dtauW(rknotsSpoly))
knotsS = 0.0_rk

CALL LINSPACE(Sbound(1), Sbound(2), rknotsS, knotsS)
CALL LOGSPACE(knotsS(2)/logbound, Sbound(2), rknotsS - 1, knotsS(2:rknotsS))

CALL SPLHS(knotsS, rknotsSint, indicator, LW, UW, dtauW)

ALLOCATE(csV(4, rknotsSpoly), csE1(4, rknotsSpoly))

csV = 0.0_rk; csE1 = 0.0_rk

ALLOCATE(V(rknotsS), E1(rknotsS))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!										!
!	Equilibrium Price Determination		!
!										!			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Bisect for p = Marginal Utility of household consumption in stationary equilibrium
ALLOCATE(SJ(Jmax), AlphaJ(Jmax), ThetaJ(Jmax), MYJ(Jmax), E1J(Jmax), ZIJraw(Jmax), &
        NYJ(Jmax), CYJ(Jmax), ECostJ(Jmax), RWJ(Jmax), samebins(Jmax))

SJ = 0.0_rk; AlphaJ = SJ; ThetaJ = SJ; MYJ = SJ; NYJ = SJ; ECostJ = SJ

distance = 1.0_rk

! Price bounds for steady state equilibrium

DO 

WRITE(*, '(1X,A, F9.6,A)', ADVANCE = 'NO') ' pl suggested = ', plsuggest,' enter p lower bound: '
READ(*,*) plow

WRITE(*, '(1X,A, F9.6,A)', ADVANCE = 'NO') ' ph suggested = ', phsuggest,' enter p upper bound: '
READ(*,*) phigh

CALL CPU_TIME(begintime)

WRITE(*,*) ' '

! Determine stationary equilibrium if marginal utility of households' consumption is plow
CALL EconomyQ(alpha, thetam, thetan, Atech, delta, Apref, beta, &
			  storage, precisionv, precisiond, precisionc, rknotsS, rknotsSpoly, &
			  rknotsSint, knotsS, plow, zbar, LW, UW, dtauW, csV, csE1, SVECV, SVECE1, &
			  indicator, Jmax, abeta, bbeta, zibound, masspoint, outputcost, storeq, &
			  storeh, distance, V, E1, kstar, Nkstar, Sstar, M, Csupply, maxJ, SJ, &
			  AlphaJ, ThetaJ, MYJ, NYJ, CYJ, q, E1J, ZIJraw, ECostJ, Estar, RWJ, samebins)

Tplow = 1.0_rk/Csupply

! Determine stationary equilibrium if marginal utility of households' consumption is phigh
CALL EconomyQ(alpha, thetam, thetan, Atech, delta, Apref, beta, &
			  storage, precisionv, precisiond, precisionc, rknotsS, rknotsSpoly, &
			  rknotsSint, knotsS, phigh, zbar, LW, UW, dtauW, csV, csE1, SVECV, SVECE1, &
			  indicator, Jmax, abeta, bbeta, zibound, masspoint, outputcost, storeq, &
			  storeh, distance, V, E1, kstar, Nkstar, Sstar, M, Csupply, maxJ, SJ, &
			  AlphaJ, ThetaJ, MYJ, NYJ, CYJ, q, E1J, ZIJraw, ECostJ, Estar, RWJ, samebins)

Tphigh = 1.0_rk/Csupply

IF ((plow-Tplow)*(phigh-Tphigh).GT.0.0_rk) THEN
	WRITE(*,*) ' price bounds do not bisect a steady state equilibrium '
	WRITE(*,'(1X,A,F10.4, A, F13.4)', ADVANCE = 'YES') ' plow =  ', plow,  ' T(plow) or 1/C(plow) =  ', Tplow
	WRITE(*,'(1X,A,F10.4, A, F13.4)', ADVANCE = 'YES') ' phigh = ', phigh, ' T(phigh) or 1/C(phigh) = ', Tphigh
ELSE 
	s1 = 0_ik
	EXIT
END IF

END DO


DO

distance = phigh - plow
IF (distance.LT.precisionc*100) THEN
	EXIT
ELSE
	s1 = s1 + 1_ik
END IF

p = (plow + phigh)/2.0_rk

CALL EconomyQ(alpha, thetam, thetan, Atech, delta, Apref, beta, &
			  storage, precisionv, precisiond, precisionc, rknotsS, rknotsSpoly, &
			  rknotsSint, knotsS, p, zbar, LW, UW, dtauW, csV, csE1, SVECV, SVECE1, &
			  indicator, Jmax, abeta, bbeta, zibound, masspoint, outputcost, storeq, &
			  storeh, distance, V, E1, kstar, Nkstar, Sstar, M, Csupply, maxJ, SJ, &
			  AlphaJ, ThetaJ, MYJ, NYJ, CYJ, q, E1J, ZIJraw, ECostJ, Estar, RWJ, samebins)

Tp0 = 1.0_rk/Csupply

IF ((p - Tp0).LT.0.0_rk) THEN
	plow = p
ELSE
	phigh = p
END IF

END DO

! Is the storage cost a fraction of the value of inventories, q*Sf?
IF (storeq.EQ.'y') THEN
	storage1 = storage*q
ELSE
	storage1 = storage
END IF

GDP = delta*kstar + Csupply
realrate = 1.0_rk/beta - 1.0_rk
cyratio = Csupply/GDP; kyratio = kstar/GDP; myratio = q*M/GDP

! total adjustment costs
adjustors = DOT_PRODUCT(AlphaJ, ThetaJ)

! Adjusted for output cost possibility and proper introduction of ECostJ
Nm = DOT_PRODUCT(NYJ(1:maxJ+1),(/adjustors, ThetaJ(1:maxJ)/))
AdCost = DOT_PRODUCT(ThetaJ(1:maxJ), ECostJ(1:maxJ))
Totalhours = Nkstar + (1.0_rk - outputcost)*AdCost + Nm

CALL CPU_TIME(endingtime); endingtime = endingtime - begintime

WRITE(*,*) ' ' 
WRITE(*,'(1X,4(A,F5.3))', ADVANCE = 'NO') ' alpha =  ', alpha, ' thetam = ', thetam, ' thetan = ', thetan
WRITE(*,'(1X,4(A,F5.3))', ADVANCE = 'YES')  ' abeta = ', abeta, ' bbeta = ', bbeta, ' zibound = ', zibound

WRITE(*,'(1X, 4(A,F5.3), A, F8.4, A, A)', ADVANCE = 'YES')  ' zbar = ', zbar, &
' delta = ', delta, ' beta = ', beta, ' Apref = ', Apref, ' storage = ', storage, ' storeh = ', storeh

WRITE(*,*) '  ' 

WRITE(*, FMT = '(1X, A, F9.5, A, F9 .5, A, F6.4, A, F6.4,A, F6.4, A, A)', ADVANCE = 'YES') &
' p0 = ', p, ' T(p0) = ', 1.0_rk/Csupply, ' (q, M) = (', q,',', M,') total adj. cost: ', AdCost, ' file: ', datafile

WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*, FMT = '(1X, A)', ADVANCE = 'NO') '    S*,Sj: '
WRITE(*, FMT = '(1X, F9.5)', ADVANCE = 'NO') Sstar
DO si = 1, maxJ; WRITE(*, FMT = '(1X, F9.5)', ADVANCE = 'NO') SJ(si); END DO

WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '
!WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*, FMT = '(1X, A)', ADVANCE = 'NO') '       Mj: '
DO si = 1, maxJ+1; WRITE(*, FMT = '(1X,F9.5)', ADVANCE = 'NO') MYJ(si); END DO
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*, FMT = '(1X, A)', ADVANCE = 'NO') '       Nj: '
DO si = 1, maxJ+1; WRITE(*, FMT = '(1X,F9.5)', ADVANCE = 'NO') NYJ(si); END DO
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*, FMT = '(1X, A)', ADVANCE = 'NO') '      CYj: '
DO si = 1, maxJ+1; WRITE(*, FMT = '(1X,F9.5)', ADVANCE = 'NO') CYJ(si); END DO
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*, FMT = '(1X, A)', ADVANCE = 'NO') '      RWj: '
DO si = 1, maxJ+1; WRITE(*, FMT = '(1X,F9.5)', ADVANCE = 'NO') RWJ(si); END DO
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

!WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*, FMT = '(1X, A)', ADVANCE = 'NO') '   alphaj:           '
DO si = 1, maxJ; WRITE(*, FMT = '(1X,F9.5)', ADVANCE = 'NO') AlphaJ(si); END DO
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*, FMT = '(1X, A)', ADVANCE = 'NO') '   thetaj:           '
DO si = 1, maxJ; WRITE(*, FMT = '(1X,F9.5)', ADVANCE = 'NO') ThetaJ(si); END DO
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*, FMT = '(1X, A)', ADVANCE = 'NO') ' prodshare '
WRITE(*, FMT = '(1X, F9.5)', ADVANCE = 'NO') adjustors
DO si = 1, maxJ; WRITE(*, FMT = '(1X,F9.5)', ADVANCE = 'NO') (1.0_rk - AlphaJ(si))*ThetaJ(si); END DO
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

!WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*, FMT = '(1X, A)', ADVANCE = 'NO') '   Ea,E1j: '   
WRITE(*, FMT = '(1X, F9.5)', ADVANCE = 'NO') Estar 
DO si = 1, maxJ; WRITE(*, FMT = '(1X,F9.5)', ADVANCE = 'NO') E1J(si); END DO
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*, FMT = '(1X, A)', ADVANCE = 'NO') ' zitildej:           '
DO si = 1, maxJ; WRITE(*, FMT = '(1X,F9.5)', ADVANCE = 'NO') ZIJraw(si); END DO
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*, FMT = '(1X, A)', ADVANCE = 'NO') '   s = 0 :          '
DO si = 1, maxJ; WRITE(*, FMT = '(1X,I8)', ADVANCE = 'NO') samebins(si); END DO
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*, FMT = '(1X, A)', ADVANCE = 'NO') '   ECostj:           '
DO si = 1, maxJ; WRITE(*, FMT = '(1X,F9.5)', ADVANCE = 'NO') ECostJ(si); END DO
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*,'(1X, 7(A, F9.5))', ADVANCE = 'YES') ' Nk: ', Nkstar, ' Nm: ', Nm, &
' hours: ', Totalhours, &
 ' q.r. rate: ', realrate, ' I/K: ', delta, ' storage1 = ', storage1
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '
WRITE(*,'(1X,A, F8.4, A, F8.4, A, F8.4, A, F8.4)', ADVANCE = 'NO') ' K: ', kstar, &
' SM: ', DOT_PRODUCT(ThetaJ,SJ), ' C: ', Csupply, ' GDP: ', GDP
WRITE(*,'(1X,3(A,F8.4))', ADVANCE = 'YES')  ' C/Y: ', cyratio, ' K/Y: ', kyratio, ' QM/Y: ', myratio
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '
WRITE(*,'(1X, 2(A, F8.4))', ADVANCE = 'YES') ' K*SM = ', kstar*DOT_PRODUCT(ThetaJ,SJ), '  Inventory/Sales Ratio = Q*SM/GDP = ',  q*DOT_PRODUCT(ThetaJ,SJ)/GDP
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(*,'(1X, (A, F8.4))', ADVANCE = 'YES') '  C + delta*k = ', Csupply + delta*kstar
	
WRITE(*, FMT = '(1X,A)', ADVANCE = 'YES') ' '	

! OUTPUT TO RESULT FILE
OPEN(UNIT = 30, FILE = resultfile, ACTION = 'READWRITE', STATUS = 'REPLACE')	

WRITE(30, FMT = '(A)', ADVANCE = "YES") "  "

WRITE(30,*) ' ' 
WRITE(30,'(1X,3(A,F5.3))', ADVANCE = 'NO') ' alpha =  ', alpha, ' thetam = ', thetam, ' thetan = ', thetan
WRITE(30,'(1X,4(A,F5.3))', ADVANCE = 'YES')  ' abeta = ', abeta, ' bbeta = ', bbeta, ' zibound = ', zibound, ' p0 = ', masspoint

WRITE(30,'(1X, 4(A,F5.3), A, F8.4, A, A)', ADVANCE = 'YES')  ' zbar = ', zbar, &
' delta = ', delta, ' beta = ', beta, ' Apref = ', Apref, ' storage = ', storage, ' storeh = ', storeh

WRITE(30,*) '  ' 

WRITE(30, FMT = '(1X, A, F9.5, A, F9 .5, A, F6.4, A, F6.4,A, F6.4, A, A)', ADVANCE = 'YES') &
' p0 = ', p, ' T(p0) = ', 1.0_rk/Csupply, ' (q, M) = (', q,',', M,') total adj. cost: ', AdCost, ' file: ', datafile

WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X, A)', ADVANCE = 'NO') '    S*,Sj: '
WRITE(30, FMT = '(1X, F8.4)', ADVANCE = 'NO') Sstar
DO si = 1, maxJ; WRITE(30, FMT = '(1X, F8.4)', ADVANCE = 'NO') SJ(si); END DO

WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '
!WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X, A)', ADVANCE = 'NO') '       Mj: '
DO si = 1, maxJ+1; WRITE(30, FMT = '(1X,F8.4)', ADVANCE = 'NO') MYJ(si); END DO
WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X, A)', ADVANCE = 'NO') '       Nj: '
DO si = 1, maxJ+1; WRITE(30, FMT = '(1X,F8.4)', ADVANCE = 'NO') NYJ(si); END DO
WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X, A)', ADVANCE = 'NO') '      CYj: '
DO si = 1, maxJ+1; WRITE(30, FMT = '(1X,F8.4)', ADVANCE = 'NO') CYJ(si); END DO
WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

!WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X, A)', ADVANCE = 'NO') '   alphaj:          '
DO si = 1, maxJ; WRITE(30, FMT = '(1X,F8.4)', ADVANCE = 'NO') AlphaJ(si); END DO
WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X, A)', ADVANCE = 'NO') '   thetaj:          '
DO si = 1, maxJ; WRITE(30, FMT = '(1X,F8.4)', ADVANCE = 'NO') ThetaJ(si); END DO
WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X, A)', ADVANCE = 'NO') ' prodshare '
WRITE(30, FMT = '(1X, F8.4)', ADVANCE = 'NO') adjustors
DO si = 1, maxJ; WRITE(30, FMT = '(1X,F8.4)', ADVANCE = 'NO') (1.0_rk - AlphaJ(si))*ThetaJ(si); END DO
WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

!WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X, A)', ADVANCE = 'NO') '      E1j:          '
DO si = 1, maxJ; WRITE(30, FMT = '(1X,F8.4)', ADVANCE = 'NO') E1J(si); END DO
WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X, A)', ADVANCE = 'NO') ' zitildej:          '
DO si = 1, maxJ; WRITE(30, FMT = '(1X,F8.4)', ADVANCE = 'NO') ZIJraw(si); END DO
WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X, A)', ADVANCE = 'NO') '   s = 0 :          '
DO si = 1, maxJ; WRITE(30, FMT = '(1X,I8)', ADVANCE = 'NO') samebins(si); END DO
WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X, A)', ADVANCE = 'NO') '   ECostj:          '
DO si = 1, maxJ; WRITE(30, FMT = '(1X,F8.4)', ADVANCE = 'NO') ECostJ(si); END DO
WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X,A, F12.8)', ADVANCE = 'YES') ' Sum of Theta = ', SUM(ThetaJ)

WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30,'(1X, 4(A, F8.4))', ADVANCE = 'YES') ' Nk: ', Nkstar, ' Nm: ', Nm, &
' hours: ', Totalhours, ' storage1 = ', storage1
WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '
WRITE(30,'(1X,A, F8.4, A, F8.4, A, F8.4, A, F8.4)', ADVANCE = 'NO') ' K: ', kstar, &
' SM: ', DOT_PRODUCT(ThetaJ,SJ), ' C: ', Csupply, ' GDP: ', GDP
WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '
WRITE(30,'(1X, 1(A, F8.4))', ADVANCE = 'YES') ' Inventory/Sales Ratio = Q*SM/GDP = ',  q*DOT_PRODUCT(ThetaJ,SJ)/GDP
WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '

WRITE(30, FMT = '(1X,A)', ADVANCE = 'YES') ' '	


CLOSE(30)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!																		!
!	Record Initial Condition for Dynamic Inventory Model's Outer Loop	!
!																		!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(UNIT = 57, FILE = startfile, ACTION = 'WRITE', STATUS = 'REPLACE')

WRITE(57,*) Jmax

DO si = 1, Jmax
	WRITE(57,*) SJ(si)
END DO

DO si = 1, Jmax
	WRITE(57,*) ThetaJ(si)
END DO

WRITE(57,*) kstar

WRITE(57,*) indicator

DO si = 1, 2
WRITE(57,*) SVECV(si,1)
END DO

DO si = 1,2
WRITE(57,*) SVECE1(si,1)
END DO

WRITE(57,*) logbound

CLOSE(57)

WRITE(*, FMT = '(1X,A,F8.2, A, A)', ADVANCE = 'YES') ' elapsed time: ', endingtime, ' steady state saved as initial condition in ', startfile


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!	Record Data for the Impulse Response Programs		!
!														!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

resultfile(1:30) = '                    ';
resultfile(1:7) = 'trans'
si = LEN_TRIM(startfile)
resultfile(6:5+si) = startfile(1:si)
resultfile(6+si:9+si) = '.txt'

OPEN(UNIT = 57, FILE = resultfile, ACTION = 'WRITE', STATUS = 'REPLACE')

! Parameters of the model
WRITE(57,*) ' alpha    =		', alpha
WRITE(57,*) ' thetam   =		', thetam
WRITE(57,*) ' thetan   =		', thetan
WRITE(57,*) ' Atech    =		', Atech
WRITE(57,*) ' delta    =		', delta
WRITE(57,*) ' Apref    =		', Apref
WRITE(57,*) ' beta     = 		', beta
WRITE(57,*) ' abeta    = 		', abeta
WRITE(57,*) ' bbeta    = 		', bbeta
WRITE(57,*) ' zibound    = 		', zibound
WRITE(57,*) ' p0       =        ', masspoint

WRITE(57,*) ' '
WRITE(57,*) ' outputcost  =		', outputcost
WRITE(57,*) ' storage	  =		', storage
WRITE(57,*) ' storage1	  =		', storage1
WRITE(57,*) ' storeh	  =		', storeh
WRITE(57,*) ' precisions  =		', precisiond
WRITE(57,*) ' precisionm  =		', precisionc
WRITE(57,*) ' Jmax		  =		', Jmax
WRITE(57,*) ' NS   		  =		', NS
WRITE(57,*) '  '

WRITE(57,*) ' Steady State Values '

WRITE(57,*) ' pss	   =		', p
WRITE(57,*) ' Qss      =		', p*q
WRITE(57,*) ' kss	   =		', kstar

WRITE(57,*) '  '
WRITE(57,*) ' knots on S and V(S)  '

WRITE(57, FMT = '(1X, A)', ADVANCE = 'NO') '(/ '
DO si = 1, NS; WRITE(57, FMT = '(1X, F12.8, A)', ADVANCE = 'NO') knotsS(si), ', '; END DO
WRITE(57, FMT = '(1X, A)', ADVANCE = 'NO') '/) '
WRITE(57,*) '  '

WRITE(57, FMT = '(1X, A)', ADVANCE = 'NO') '(/ '
DO si = 1, NS; WRITE(57, FMT = '(1X,F12.8, A)', ADVANCE = 'NO') V(si), ', '; END DO
WRITE(57, FMT = '(1X, A)', ADVANCE = 'NO') '/) '
WRITE(57,*) '  '

WRITE(57,*) '  '
WRITE(57,*) ' Steady State Svec and MUvec  '

WRITE(57, FMT = '(1X, A)', ADVANCE = 'NO') '(/ '
DO si = 1, maxJ; WRITE(57, FMT = '(1X, F12.8, A)', ADVANCE = 'NO') SJ(si), ', '; END DO
WRITE(57, FMT = '(1X, A)', ADVANCE = 'NO') '/) '
WRITE(57,*) '  '

WRITE(57, FMT = '(1X, A)', ADVANCE = 'NO') '(/ '
DO si = 1, maxJ; WRITE(57, FMT = '(1X,F12.8, A)', ADVANCE = 'NO') ThetaJ(si), ', '; END DO
WRITE(57, FMT = '(1X, A)', ADVANCE = 'NO') '/) '
WRITE(57,*) '  '

CLOSE (57)



!! ************************************************************************ !!
!!																			!!
!!							INTERNAL SUBROUTINES							!!
!!																			!!
!! ************************************************************************ !!

CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!							!!
!!							!!	
!!							!!	
!!  	Economy (p,q)		!!  	 
!!							!!	
!!							!!		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Essential subroutine to determine steady state behavior of inventory firms
! given p and q.
SUBROUTINE EconomyQ(alpha, thetam, thetan, Atech, delta, Apref, beta, &
					storage0, precisionv, precisiond, precisionc, rknotsS, rknotsSpoly, &
					rknotsSint, knotsS, p, zbar, LW, UW, dtauW, csV, csE1, SVECV, SVECE1, &
					indicator, Jmax, abeta, bbeta, zibound, masspoint, outputcost, storeq, &
					storeh, distancep, V, E1, kstar, Nkstar, Sstar, X, Ctotal, maxJ, SJ, &
					AlphaJ, ThetaJ, MYJ, NYJ, CYJ, q, E1S, ZIJraw, ECostJ, Estar, RWJ, samebins)


IMPLICIT NONE

INTEGER(ik):: rknotsS, rknotsSint, rknotsSpoly, Jmax, maxJ, si, s1, gooble, &
			  samebins(Jmax), maxJs, gooble2
REAL(rk):: alpha, thetam, thetan, Atech, delta, Apref, beta, precisionv, &
		   precisiond, precisionc, knotsS(rknotsS), q, p, zbar, LW(1, rknotsSint-1), &
		   UW(2, rknotsSint), dtauW(rknotsSpoly), csV(4, rknotsSpoly), Sval, costterm, &
		   csE1(4, rknotsSpoly), SVECV(2,1), SVECE1(2,1), abeta, bbeta, zibound, masspoint, &
		   V(rknotsS), E1(rknotsS), kstar, Nkstar, Sstar, expterm, gz, chi, &
		   dv, de, distance, E1S(Jmax), X, Ea, msol, nsol, cysol, &
		   zitilde, Vf(1,1), Sf(1), Ctotal, alphai, SJ(Jmax), AlphaJ(Jmax), ThetaJ(Jmax), &
		   CYJ(Jmax), outputcost, MYJ(Jmax), ZIJ(Jmax), adjustors, distancep, ZIJraw(Jmax), &
		   NYJ(Jmax), ECostJ(Jmax), storage, Estar, RWJ(Jmax), storage0, storageh
REAL(rk), ALLOCATABLE:: TV(:), TE1(:)
CHARACTER(20):: indicator, storeq, storeh
INTENT(IN):: rknotsS, alpha, thetam, thetan, Atech, delta, Apref, beta, &
			 precisionv, precisiond, precisionc, knotsS, p, zbar, LW, storage0, &
		     UW, dtauW, Jmax, indicator, distancep, SVECV, SVECE1, storeq, storeh, &
			 abeta, bbeta, zibound, masspoint
INTENT(OUT):: V, E1, kstar, Nkstar, Sstar, X, maxJ, SJ, AlphaJ, ThetaJ, MYJ, q, E1S, &
			  ZIJraw, Ctotal, ECostJ, Estar, RWJ
INTENT(INOUT):: csV, csE1

! Initial value functions
V = 10.0_rk
E1 = 0.0_rk

q = ((1.0_rk - beta*(1.0_rk - delta))/(beta*alpha))**alpha
q = q*((Apref/(1.0_rk - alpha))**(1.0_rk - alpha))
q = (q/zbar)*(p**(alpha - 1.0_rk))

IF (storeq.EQ.'y') THEN
	storage = storage0*q
ELSE
	storage = storage0
END IF

IF(storeh.EQ.'y') THEN	! Return storage cost to households
	storageh = 0.0_rk
ELSE
	storageh = 1.0_rk
END IF

! numerical terms
gz = Atech*zbar 

expterm = thetam/(1.0_rk - thetan)				
chi = (thetan/Apref)**thetan
chi = (gz*chi)**(1.0_rk/(1.0_rk - thetan))
chi = (1.0_rk - thetan)*chi

DO si = 1, NS							! initial value function
	E1(si) = p**(1.0_rk/(1.0_rk - thetan))
	E1(si) = E1(si)*chi*(knotsS(si)**expterm)
END DO

ALLOCATE(TV(rknotsS), TE1(rknotsS))

V = 0.0_rk
distance = 2*precisionv
TV =  0.0_rk
TE1 = 0.0_rk


!!!!!!!!!!!!!!!!!!!!!!!!!
!						!
!	Iterate on E1 and V	!
!						!
!!!!!!!!!!!!!!!!!!!!!!!!!

s1 = 0_ik

DO

	IF (distance.LT.precisionv) THEN
		EXIT
	ELSE
		s1 = s1 + 1_ik
	END IF

ZIJ = zibound; SJ = 0.0_rk; CYJ = 0.0_rk; MYJ = 0.0_rk; E1S = 0.0_rk
ZIJraw = 0.0_rk; NYJ = 0.0_rk

! Fit a spline to update polynomial approximations to E1 and V

CALL SPpp(E1, rknotsSint, 1, LW, UW, dtauW, indicator, SVECE1, csE1(:,:))
CALL SPpp(V, rknotsSint, 1, LW, UW, dtauW, indicator, SVECV, csV(:,:))

! Locate S*, the choice of target inventory stock by adjusting firms using E1(S)

CALL InvSstar(rknotsS, p*q, precisionc, knotsS, csE1, (/knotsS(1), knotsS(rknotsS)/), Sstar, Ea)

! Now determine thresholds using firms alternative valuation of V(S').  Those 
! firms that adjusted last period, will have the largest beginning of this period 
! S, and this will be S(1) = S* - M(S*).  We need to compute this S(1) by 
! determining M(S*).  No firm begins the period with S*
CALL InvM(rknotsS, thetam, thetan, Atech, p, beta, precisiond, Apref, storage, &
			 storageh, knotsS, csV, (/0.0_rk, Sstar/), zbar, Sstar,  msol, nsol, cysol, E1S(1))

! Firms that adjust to S* will produce cysol and use msol.
CYJ(1) = cysol
MYJ(1) = msol
NYJ(1) = nsol

! However, the beginning of the period distribution of S will start with SJ(1).
Estar = Ea
SJ(1) = Sstar - msol
RWJ(1) = p*cysol - Apref*nsol - (1.0_rk - storageh)*p*storage*SJ(1)

! Note largest feasible history is set to Jmax - 1
maxJ = Jmax - 1_ik

gooble = 0_ik; gooble2 = 0_ik; samebins = 0_ik

! Labor or output cost
costterm = (1.0_rk - outputcost)*Apref + outputcost*p

DO si = 1, Jmax
	
	IF (SJ(si).GT.precisionc) THEN

	CALL InvM(rknotsS, thetam, thetan, Atech, p, beta, precisiond, Apref, storage, &
				 storageh, knotsS, csV, (/0.0_rk, SJ(si)/), zbar, SJ(si),  msol, nsol, cysol, E1S(si))
		gooble2 = 0_ik

	ELSE
		
		msol = SJ(si); Sf(1) = 0.0_rk; CALL SPeval(csV, knotsS, rknotsS - 2_ik, 1, Sf, 1, 0, Vf)
		nsol = thetan*p*gz*(msol**thetam)/Apref
		nsol = nsol**(1.0_rk/(1.0_rk - thetan))
		cysol = gz*(msol**thetam)*(nsol**thetan)
		E1S(si) = p*cysol - Apref*nsol + beta*Vf(1,1)
		gooble2 = 1_ik

	END IF
	
	! This is equation (23) on page 6.					
	zitilde = (p*q*SJ(si) - E1S(si) + Ea)/costterm		 
	ZIJraw(si) = zitilde
	
	IF (si.LT.Jmax) THEN

		IF (gooble2.NE.1_ik) THEN
			SJ(si+1) = SJ(si) - msol
		ELSE
			SJ(si+1) = 0.0_rk
			samebins(si+1) = 1_ik			
		END IF

		CYJ(si+1) = cysol
		MYJ(si+1) = msol
		NYJ(si+1) = nsol

		RWJ(si+1) = p*cysol - Apref*nsol - (1.0_rk - storageh)*p*storage*SJ(si+1)

		ZIJ(si) = DMIN1(DMAX1(0.0_rk, zitilde), zibound)
	ELSE
		ZIJ(si) = zibound
	END IF

	IF (zitilde.GE.zibound.AND.gooble.NE.1_ik) THEN	! All firms have adjusted by maxJ
		maxJ = si
		gooble = 1_ik
		!EXIT											
	END IF
			
END DO

! maxJ .LE. Jmax - 1, thusfar it is .EQ. unless zitile.GE.Bc
maxJs = Jmax + 1 - SUM(samebins)
maxJ = MIN(maxJ, maxJs)

! Now update E1 and V
DO si = 1, NS
	Sval = knotsS(si)


	IF (Sval.GT.precisionc) THEN

		CALL InvM(rknotsS, thetam, thetan, Atech, p, beta, precisiond, Apref, storage, &
					 storageh, knotsS, csV, (/0.0_rk, Sval/), zbar, Sval,  msol, nsol, cysol, TE1(si))
	ELSE
		
		msol = Sval; Sf(1) = 0.0_rk; CALL SPeval(csV, knotsS, rknotsS - 2_ik, 1, Sf, 1, 0, Vf)
		nsol = thetan*p*gz*(msol**thetam)/Apref
		nsol = nsol**(1.0_rk/(1.0_rk - thetan))
		cysol = gz*(msol**thetam)*(nsol**thetan)
		TE1(si) = p*cysol - Apref*nsol + beta*Vf(1,1)
		
	END IF


	! adjusted for labor or output cost
	zitilde = (p*q*Sval - TE1(si) + Ea)/costterm
	zitilde = DMIN1(DMAX1(0.0_rk, zitilde), zibound)

	! Including any integer as the 5th argument returns the probability of zi.le.zitilde
	alphai = CostFunction(zitilde, abeta, bbeta, zibound, masspoint, 0_ik)

	TV(si) = alphai*(p*q*Sval + Ea) - costterm*CostFunction(zitilde, abeta, bbeta, zibound, masspoint)
	TV(si) = TV(si) + (1.0_rk - alphai)*TE1(si)

END DO	

de = MAXVAL(DABS(TE1-E1)); dv = MAXVAL(DABS(TV - V))
distance = DMAX1(de,dv)

IF (s1.GT.1750) THEN
	WRITE(*, '(A, 1X, I6, A, E8.2, A, F6.3, A, E8.2, A, E8.2, A, I2, A, F8.4, A, F8.4)', ADVANCE = 'YES') &
	' RAISE S* (', s1, ') (gap, p): (', distancep, ', ', p, ')  |E1| = ', de, ' |V| = ', dv, ' J: ', maxJ, '   S*: ', Sstar, '   alphai: ', alphai

END IF

V = TV
E1 = TE1

END DO ! end of contraction loop on E1 and V.

WRITE(*, '(1X, I6, A, E8.2, A, 2(F8.4), A, E12.4, A, I4, A, F8.4, A, E10.2)', ADVANCE = 'YES') &
s1, ') (gap, p, q): (', distancep, ', ', p, q, ')  |E1/V|: ', distance, ' J: ', maxJ, '  S*: ', Sstar, '  E1a(S*): ', Ea

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							!
!	The MU distribution		!
!							!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

AlphaJ = 0.0_rk; ThetaJ = 0.0_rk; 

! Now determine distribution mu, first determine adjustment rates, alpha(j).
! Remember that adjustment occurs at the beginning of each period, before 
! production, hence alpha(1) = 0 as firms here are exactly those that 
! have adjusted to S*.
DO si = 1, Jmax
	AlphaJ(si)  = CostFunction(ZIJ(si), abeta, bbeta, zibound, masspoint, 0_ik)
	ECostJ(si) = CostFunction(ZIJ(si), abeta, bbeta, zibound, masspoint)
END DO
					
ThetaJ(1) = 1.0_rk

DO si = 1, maxJ - 1	! Note maxJ is defined to be no more than Jmax - 1
	ThetaJ(si+1) = (1.0_rk - AlphaJ(si))*ThetaJ(si)
END DO

! Note this adjustment is consistent whether AlphaJ(maxJ) = 1 or not.
ThetaJ(maxJ) = ThetaJ(maxJ)/AlphaJ(maxJ)

Xnumadjust = 1.0_rk/SUM(ThetaJ)
ThetaJ = 0.0_rk
ThetaJ(1) = Xnumadjust

DO si = 1, maxJ - 1
	ThetaJ(si+1) = (1.0_rk - AlphaJ(si))*ThetaJ(si)
END DO

ThetaJ(maxJ) = ThetaJ(maxJ)/AlphaJ(maxJ)

X = SUM((Sstar - SJ(1:maxJ))*AlphaJ(1:maxJ)*ThetaJ(1:maxJ))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!												!
!	Determining K and Nk so that X = F(K,Nk)	!
!												!		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Determine steady state capital K.  Given we use above X, which is Xd, this implies 
! that Xd = Xs and the levels of K, and later Nk, will yield production of X above.  

kstar = (p*q*(1.0_rk - alpha)/Apref)**((1.0_rk - alpha)/alpha)
kstar = (zbar**(1.0_rk/alpha))*kstar
kstar = X/kstar

! Given kstar and q, Nk is given by equation 6. 
Nkstar = q*zbar*(1.0_rk - alpha)
Nkstar = Nkstar/(Apref/p)
Nkstar = Nkstar**(1.0_rk/alpha)
Nkstar = Nkstar*kstar

! Adjusting plants reach S* before production and produce CYJ(1).  
! Plants with S(j) who do not adjust produce CYJ(j+1).  Must adjust 
! Ctotal for output costs if necessary.  Moreover, note that below
! both consumption and investment are final goods and that CYJ is 
! thus the vector of output produced by final goods firms and not
! the vector of consumption.
adjustors = SUM(AlphaJ*ThetaJ)

Ctotal = SUM((/adjustors, (1.0_rk - AlphaJ(1:maxJ))*ThetaJ(1:maxJ)/)*CYJ(1:maxJ+1)) - delta*kstar - &
		 outputcost*DOT_PRODUCT(ThetaJ(1:maxJ), ECostJ)

END SUBROUTINE EconomyQ


!!!!!!!!!!!!!!
!!!!!!!!!!!!!!
!!			!!
!!	SETUP	!!
!!			!!
!!!!!!!!!!!!!!
!!!!!!!!!!!!!!

SUBROUTINE SETUP(alpha, thetam, thetan, Atech, delta, abeta, bbeta, zibound, masspoint, &
				 outputcost, zbar, Apref, beta, storage, precisionv, precisiond, precisionc, NS, &
				 Sbound, Jmax, datafile, resultfile, startfile, pl0, ph0, logbound, indicator, &
				 SVECV, SVECE1)  

INTEGER(ik):: NS, i, what, Jmax
REAL(rk):: alpha, thetam, thetan, Atech, delta, beta, storage, &
		   precisionv, precisiond, precisionc, Apref, switch, Sbound(2), abeta, bbeta, zibound, masspoint, &
		   zbar, pl0, ph0, logbound, outputcost, SVECV(2,1), SVECE1(2,1)
CHARACTER(30):: datafile, resultfile, startfile, indicator

INTENT(OUT):: alpha, thetam, thetan, Atech, delta, beta, outputcost, precisionc, &
			  precisionv, precisiond, datafile, resultfile, Sbound, Apref, NS, pl0, ph0, indicator, &
			  SVECV, SVECE1, abeta, bbeta, zibound, masspoint


DO 

	WRITE(*,FMT = '(1X,A)', ADVANCE = 'NO') ' enter parameter file: '
	READ(*,*) datafile
	WRITE(*,*) ' '
	
	OPEN (UNIT = 57, FILE = datafile, STATUS = 'OLD', ACTION = 'READ', IOSTAT = what, ERR = 100)
	
	! The subroutine goes to line 100, right below, if there's a file opening error.  
	100	IF (what.NE.0_ik) THEN 
			WRITE(*,'(1X,A,A)', ADVANCE = 'YES') ' Hey, I think I could not locate ', datafile
			WRITE(*,*) ' I got a file I/O error ', what
			what = 0_ik
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
READ(57,*) zbar											! mean level of z
READ(57,*)
READ(57,*) abeta										! beta distribution exponent 1 parameter
READ(57,*) 
READ(57,*) bbeta										! beta distribution exponent 2 parameter
READ(57,*)
READ(57,*) zibound										! cost distribution upper bound
READ(57,*)
READ(57,*) masspoint									! probability of zero cost
READ(57,*)
READ(57,*) outputcost									! Are adjustment costs in units of output? 1 for yes
READ(57,*)
READ(57,*) beta											! discount factor
READ(57,*)
READ(57,*) Apref										! preference term for leisure
READ(57,*)
READ(57,*) storage										! storage cost for inventories of the intermediate input
READ(57,*)
READ(57,*) (Sbound(i), i = 1,2)							! knot endpoints on capital
READ(57,*)
READ(57,*)  NS											! total number of knots, including endpoints, on S	
READ(57,*)
READ(57,*)  switch										! set to 1 for second order Taylor expansions
READ(57,*)  
READ(57,*)  precisionv									! precision term for value function convergence
READ(57,*)  
READ(57,*)  precisionc									! precision term for S choice
READ(57,*) 
READ(57,*)  precisiond									! precision term for m choice
READ(57,*)  
READ(57,*)  Jmax										! maximum number of groups
READ(57,*)  
READ(57,*)  indicator									! not-a-knot, complete or natural spline endpoint condition
READ(57,*)
READ(57,*) (SVECV(i,1), i = 1,2)						! endpoint condition values for V
READ(57,*)
READ(57,*) (SVECE1(i,1), i = 1,2)						! endpoint condition values for E1
READ(57,*)
READ(57,*)  resultfile									! result file name
READ(57,*)  
READ(57,*)  startfile									! file in which to record initial conditions for dynamics
READ(57,*)  
READ(57,*)  pl0											! lower bound on p in steady state
READ(57,*)  
READ(57,*)  ph0											! upper bound on p in steady state	
READ(57,*, IOSTAT = what, ERR = 300)  											
READ(57,*, IOSTAT = what, ERR = 300) logbound								

CLOSE (57)						! close datafile

! The subroutine comes to line 300 if the file ended after ph0 and logbound was not available.
300	IF (what.NE.0_ik) THEN 
			WRITE(*,'(1X,A,I4)', ADVANCE = 'YES') ' At the end of the file, I got an I/O error = ', what
			WRITE(*,*) ' Since logspacing parameter is not available, I will linspace. '
			CLOSE (57)
			WRITE(*,*) ' '
			logbound = 0.0_rk
	END IF

END SUBROUTINE SETUP


!!!!!!!!!!!!!!!!!!!!!!!!!
!						!
!	Adjustment Cost		!
!						!
!						!
!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION CostFunction(zi, abeta, bbeta, zibound, masspoint, choice)

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


!!!!!!!!!!!!!!!!!
!				!			
!	LinSpace	!
!				!
!!!!!!!!!!!!!!!!!

SUBROUTINE LinSpace(lb, ub, gridnum, X)

INTEGER(ik):: gridnum, j1
REAL(rk):: lb, ub, X(gridnum), Y(gridnum-2_ik)

INTENT(IN):: lb, ub, gridnum
INTENT(OUT):: X

! This subroutine generates a vector which contains a grid 
! of n equally spaced elements between loweround and upper-
! bound.  Note this version produces a row vector.

DO j1 = 1_ik, gridnum - 2_ik
	Y(j1) = dble(j1)
END DO				

Y = ((ub - lb)/dble(gridnum-1_ik))*Y + lb

X(1_ik) = lb
X(2_ik:gridnum-1_ik) = Y
X(gridnum) = ub

END SUBROUTINE LinSpace

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

END PROGRAM SteadyState



!!	*****************************************************************************************************	!!
!!  *****************************************************************************************************	!!
!!																											!!
!!																											!!
!!								E X T E R N A L   S U B R O U T I N E S										!!
!!																											!!
!!																											!!			
!!	*****************************************************************************************************	!!
!!  *****************************************************************************************************	!!


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
             csS, Sbound, Sval, zval, storage
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


!!!!!!!!!!!!!!!!!
!				!			
!	LinSpace	!
!				!
!!!!!!!!!!!!!!!!!

SUBROUTINE LinSpace(lb, ub, gridnum, X)

USE KindSet

INTEGER(ik):: gridnum, j1
REAL(rk):: lb, ub, X(gridnum), Y(gridnum-2_ik)

INTENT(IN):: lb, ub, gridnum
INTENT(OUT):: X

DO j1 = 1_ik, gridnum - 2_ik
	Y(j1) = dble(j1)
END DO

Y = ((ub - lb)/dble(gridnum-1_ik))*Y + lb

X(1_ik) = lb
X(2_ik:gridnum-1_ik) = Y
X(gridnum) = ub

END SUBROUTINE LinSpace