PROGRAM OutputInv

! Program to generate data table for Khan and Thomas (2007).
!
! Programs Output (HPFILTERS)
!    Module KindSet


USE KINDSET

IMPLICIT NONE

INTEGER(ik):: simlength, Jmax, t, zi, iteration
REAL(rk):: avg(15), std(15), corry(15), lambdasmooth, simnum, corrfi, corrysi, &
			outputcost, Apref, beta, highestKM, lowestKM, &
			highestSM, lowestSM, adjustmoments(2,5)
			
REAL(rk), ALLOCATABLE:: Ksim(:), Psim(:), Qsim(:), Csim(:), Ysim(:), Xsim(:), Isim(:), NKsim(:), &
						Ssim(:,:), MUsim(:,:), Yd(:,:), Md(:,:), NCd(:,:), Alphad(:,:), S1(:),&
						adjustors(:), Msim(:), NMsim(:), GDPsim(:), THsim(:), QMsim(:), &
						YHP(:), CHP(:), IHP(:), adjustHP(:), Z1HP(:), Z2HP(:), QXsim(:),&
						THHP(:), MHP(:), QHP(:), TREND(:), S(:), NII(:), NIIHP(:), S1HP(:), &
						FS(:), FSHP(:), KHP(:), TRENDY(:), IS(:), ISHP(:), AdjCost(:), XHP(:), &
						QSf(:), QSfHP(:), DS(:), QXHP(:), QMHP(:), &
						RealWage(:), RealWageHP(:), Zsim(:,:), Targetsim(:), TargetHP(:), adjustvec(:,:), rf(:), &
						S1dist(:,:), meanS1(:), varS1(:), cvS1(:), cvS1HP(:), &
						cm1(:), Gap(:,:), adjust2(:)	 

CHARACTER(30):: resultfile, dynamicfile
EXTERNAL HPFilterS


WRITE(*,'(1X,A)', ADVANCE = 'NO') ' enter inventory model result file name: '
READ(*,*) resultfile

OPEN(UNIT = 45, FILE = resultfile, ACTION = 'READ', STATUS = 'OLD')	

READ(45,*) dynamicfile
READ(45,*) iteration
READ(45,*) simlength
READ(45,*) Jmax
READ(45,*) outputcost

ALLOCATE(Ksim(simlength+1), Psim(simlength), Qsim(simlength), Csim(simlength), &
		 Xsim(simlength), Isim(simlength), NKsim(simlength), Ssim(simlength+1, Jmax), &
		 MUsim(simlength+1, Jmax), Yd(simlength, Jmax+1), Md(simlength, Jmax+1), &
		 NCd(simlength, Jmax+1), Alphad(simlength, Jmax), GDPsim(simlength), &
		 IS(simlength), AdjCost(simlength), Ysim(simlength), THsim(simlength), &
		 Zsim(2,simlength))

DO t = 1, simlength+1; READ(45,*) Ksim(t); END DO
DO t = 1, simlength; READ(45,*) Psim(t); END DO
DO t = 1, simlength; READ(45,*) Qsim(t); END DO
DO t = 1, simlength; READ(45,*) Csim(t); END DO
DO t = 1, simlength; READ(45,*) Ysim(t); END DO
DO t = 1, simlength; READ(45,*) Xsim(t); END DO
DO t = 1, simlength; READ(45,*) Isim(t); END DO
DO t = 1, simlength; READ(45,*) NKsim(t); END DO
DO t = 1, simlength+1; DO zi = 1, Jmax; READ(45,*) Ssim(t,zi); END DO; END DO
DO t = 1, simlength+1; DO zi = 1, Jmax; READ(45,*) MUsim(t,zi); END DO; END DO
DO t = 1, simlength; DO zi = 1, Jmax + 1; READ(45,*) Yd(t,zi); END DO; END DO
DO t = 1, simlength; DO zi = 1, Jmax + 1; READ(45,*) Md(t,zi); END DO; END DO
DO t = 1, simlength; DO zi = 1, Jmax + 1; READ(45,*) NCd(t,zi); END DO; END DO
DO t = 1, simlength; DO zi = 1, Jmax; READ(45,*) Alphad(t,zi); END DO; END DO
DO t = 1, simlength; READ(45,*) AdjCost(t); END DO

CLOSE(45)
 
WRITE(*,*) ' NOTE: set varnum = 1 for Matlab and Eviews output'
WRITE(*,*)

! Read the simulation data
WRITE(*,'(1X,A)', ADVANCE = 'NO') ' Does somedata with z and Apref exist? [y/n] : '
READ(*,*) resultfile

IF (resultfile.EQ.'y') THEN

	OPEN(UNIT = 35, FILE = 'somedata', ACTION = 'READ', STATUS = 'OLD')

	DO t = 1,simlength
		READ(35,*) Zsim(1,t), Zsim(2,t)
	END DO
	READ(35,*) Apref, beta

	CLOSE(35)

ELSE

	Zsim = 1.0_rk
	Apref = 0.0_rk
	beta = 1.0_rk

END IF

ALLOCATE(adjustors(simlength), Msim(simlength), NMsim(simlength), S(simlength+1), &
		 NII(simlength), FS(simlength), QSf(simlength), DS(simlength), QMsim(simlength), &
		 QXsim(simlength), S1(simlength), &
		 RealWage(simlength), Targetsim(simlength), adjustvec(simlength, Jmax), &
		 meanS1(simlength), varS1(simlength), cvS1(simlength), S1dist(simlength, Jmax), &
		 rf(simlength-1), cm1(simlength), Gap(simlength, 1:Jmax), adjust2(simlength) )
		 


! Generate relevant series
DO t = 1, simlength	
	adjustors(t) = SUM(Alphad(t,:)*MUsim(t,:))
	NMsim(t) = SUM(MUsim(t+1,1:Jmax)*NCd(t,1:Jmax))
	THsim(t) = NKsim(t) + NMsim(t) + (1.0_rk - outputcost)*AdjCost(t)
	Msim(t) = SUM(MUsim(t+1,1:Jmax)*Md(t,1:Jmax))
	FS(t) = Csim(t) + Isim(t)
	QMsim(t) = Qsim(t)*Msim(t)
	QXsim(t) = Qsim(t)*Xsim(t)
	RealWage(t) = Apref*Csim(t)	
	adjustvec(t,1:Jmax) = (Alphad(t,:)*MUsim(t,:))/adjustors(t)
END DO

S(1) = SUM(Ssim(1,:)*MUsim(1,:))

DO t = 1, simlength
	S(t+1) = SUM(Ssim(t+1,:)*MUsim(t+1,:))
	DS(t) = S(t+1) - S(t)
	NII(t) = Qsim(t)*DS(t)
	GDPsim(t) = FS(t) + NII(t)
	IS(t) = Qsim(t)*S(t+1)/FS(t)
	Targetsim(t) = Ssim(t+1, 1) + Md(t,1)
	S1dist(t,1:Jmax) = (/Targetsim(t), Ssim(t, 1:Jmax - 1)/)
	meanS1(t) = DOT_PRODUCT(S1dist(t,1:Jmax), (/adjustors(t), (1.0_rk - Alphad(t,1:Jmax-1))*MUsim(t,1:Jmax-1)/))
	varS1(t) = DOT_PRODUCT(S1dist(t,1:Jmax)**2.0_rk, (/adjustors(t), (1.0_rk - Alphad(t,1:Jmax-1))*MUsim(t,1:Jmax-1)/))
	varS1(t) = varS1(t) - meanS1(t)**2.0_rk
	cvS1(t) = varS1(t)/(meanS1(t)**2.0_rk)
	
	Gap(t, 1:Jmax) = Targetsim(t) - Ssim(t, 1:Jmax)
	adjust2(t) = DOT_PRODUCT(Gap(t, 1:Jmax), Alphad(t,1:Jmax)*MUsim(t,1:Jmax))
	adjust2(t) = adjust2(t)/DOT_PRODUCT(Gap(t, 1:Jmax), MUsim(t,1:Jmax))
	
END DO

highestKM = MAXVAL(Ksim)
lowestKM = MINVAL(Ksim)
highestSM = MAXVAL(S)
lowestSM = MINVAL(S)

! end of the period stock of inventories as in the data
QSf = Qsim*S(2:simlength+1)

! real interest rate
DO t = 1, simlength - 1
	rf(t) = Csim(t+1)/(Csim(t)*beta) - 1.0_rk
END DO

! Production time available stock of intermediate goods	
S1 = S(1:simlength) + Xsim(1:simlength)

simnum = DBLE(simlength)
avg(1) = SUM(FS)/simnum
avg(2) = SUM(GDPsim)/simnum
avg(3) = SUM(Csim)/simnum
avg(4) = SUM(Isim)/simnum
avg(5) = SUM(THsim)/simnum
avg(6) = SUM(Ksim)/(simnum+1.0_rk)
avg(7) = SUM(adjustors)/simnum
avg(8) = SUM(NII)/simnum
avg(9) = SUM(Qsim)/simnum
avg(10) = SUM(S1)/simnum
avg(11) = SUM(Msim)/(simnum)
avg(12) = SUM(IS)/simnum
avg(13) = SUM(Xsim)/simnum
avg(14) = SUM(QSf)/simnum
avg(15) = SUM(Targetsim)/simnum

adjustmoments = 0.0_rk
adjustmoments(1,1) = SUM(adjustors)/simnum
adjustmoments(2,1) = SUM(adjust2)/simnum
adjustmoments(1,2) = STDDEV(adjustors)
adjustmoments(2,2) = STDDEV(adjust2)
adjustmoments(1,3) = MAXVAL(adjustors)
adjustmoments(2,3) = MAXVAL(adjust2)
adjustmoments(1,4) = MINVAL(adjustors)
adjustmoments(2,4) = MINVAL(adjust2)
adjustmoments(1,5) = CORRELATION(adjustors, GDPsim)
adjustmoments(2,5) = CORRELATION(adjust2, GDPsim)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute Hodrick-Prescott Filtered Statistics !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

lambdasmooth = 1600.0_rk

ALLOCATE(TREND(simlength), YHP(simlength),  &
		 CHP(simlength), IHP(simlength), adjustHP(simlength), S1HP(simlength), &
		 THHP(simlength), MHP(simlength), QHP(simlength), NIIHP(simlength), &
		 FSHP(simlength), KHP(simlength+1), TRENDY(simlength), ISHP(simlength), &
		 XHP(simlength), QSfHP(simlength), QXHP(simlength), QMHP(simlength), &
		 Z1HP(simlength), Z2HP(simlength), &
		 RealWageHP(simlength), TargetHP(simlength), cvS1HP(simlength))


YHP = DLOG(GDPsim)
CALL HPFILTERS(simlength, YHP, lambdasmooth, TREND)
YHP = YHP - TREND

TREND = GDPsim
CALL HPFILTERS(simlength, TREND, lambdasmooth, TRENDY)

Z1HP = DLOG(Zsim(1,1:simlength))
CALL HPFILTERS(simlength, Z1HP, lambdasmooth, TREND)
Z1HP = Z1HP - TREND

Z2HP = DLOG(Zsim(2,1:simlength))
CALL HPFILTERS(simlength, Z2HP, lambdasmooth, TREND)
Z2HP = Z2HP - TREND

CHP = DLOG(Csim)
CALL HPFILTERS(simlength, CHP, lambdasmooth, TREND)
CHP = CHP - TREND

IHP = Isim  			! separate method of detrending investment in case I(t)<0
CALL HPFILTERS(simlength, IHP, lambdasmooth, TREND)
IHP = (IHP - TREND)/TREND

adjustHP = adjustors	! adjustment rates are not logged because they are already rates.
CALL HPFILTERS(simlength, adjustHP, lambdasmooth, TREND)
adjustHP = adjustHP - TREND

THHP = DLOG(THsim)
CALL HPFILTERS(simlength, THHP, lambdasmooth, TREND)
THHP = THHP - TREND
MHP = DLOG(Msim)
CALL HPFILTERS(simlength, MHP, lambdasmooth, TREND)
MHP = MHP - TREND
XHP = DLOG(Xsim)
CALL HPFILTERS(simlength, XHP, lambdasmooth, TREND)
XHP = XHP - TREND
QHP = DLOG(Qsim)
CALL HPFILTERS(simlength, QHP, lambdasmooth, TREND)
QHP = QHP - TREND


NIIHP = NII/GDPsim				! separate method of detrending net inventory investment
CALL HPFILTERS(simlength, NIIHP, lambdasmooth, TREND)
NIIHP = (NIIHP - TREND)

FSHP = DLOG(FS)
CALL HPFILTERS(simlength, FSHP, lambdasmooth, TREND)
FSHP = FSHP - TREND

ISHP = IS					    ! inventory to sales ratio detrended as a fraction
CALL HPFILTERS(simlength, ISHP, lambdasmooth, TREND)
ISHP = ISHP - TREND	

QSfHP = DLOG(QSf)
CALL HPFILTERS(simlength, QSfHP, lambdasmooth, TREND)
QSfHP = (QSfHP - TREND)/TREND

QXHP = Qsim*Xsim
CALL HPFILTERS(simlength, QXHP, lambdasmooth, TREND)
QXHP = (QXHP - TREND)/TRENDY

QMHP = Qsim*Msim
CALL HPFILTERS(simlength, QMHP, lambdasmooth, TREND)
QMHP = (QMHP - TREND)/TRENDY

S1HP = DLOG(S1)
CALL HPFILTERS(simlength, S1HP, lambdasmooth, TREND)
S1HP = S1HP - TREND

RealWageHP = DLOG(RealWage)		
CALL HPFILTERS(simlength, RealWageHP, lambdasmooth, TREND)
RealWageHP = RealWageHP - TREND

TargetHP = DLOG(Targetsim)
CALL HPFILTERS(simlength, TargetHP, lambdasmooth, TREND)
TargetHP = TargetHP - TREND

DEALLOCATE(TREND); ALLOCATE(TREND(simlength+1))

CALL HPFILTERS(simlength+1, DLOG(Ksim), lambdasmooth, TREND)
KHP = DLOG(Ksim) - TREND

DEALLOCATE(TREND)

std(1) = STDDEV(FSHP)	!Standard Deviations of each Series
std(2) = STDDEV(YHP)
std(3) = STDDEV(CHP)	
std(4) = STDDEV(IHP)
std(5) = STDDEV(THHP)
std(6) = STDDEV(KHP)
std(7) = STDDEV(adjustHP)	
std(8) = STDDEV(NIIHP)
std(9) = STDDEV(QHP)	
std(10) = STDDEV(S1HP)	
std(11) = STDDEV(MHP)	
std(12) = STDDEV(ISHP)
std(13) = STDDEV(XHP)
std(14) = STDDEV(QSfHP)
std(15) = STDDEV(TargetHP)	

std = std*100.0_rk

corry(1)= CORRELATION(FSHP, YHP)		  !Cross Correlations with Output
corry(2)= CORRELATION(YHP, YHP)
corry(3)= CORRELATION(CHP, YHP)
corry(4)= CORRELATION(IHP, YHP)
corry(5)= CORRELATION(THHP, YHP)
corry(6)= CORRELATION(KHP(1_ik:simlength), YHP)
corry(7)= CORRELATION(adjustHP, YHP)
corry(8)= CORRELATION(NIIHP, YHP)
corry(9)= CORRELATION(QHP, YHP)
corry(10)= CORRELATION(S1HP, YHP)
corry(11)= CORRELATION(MHP, YHP)
corry(12)= CORRELATION(ISHP, YHP)
corry(13)= CORRELATION(XHP, YHP)
corry(14)= CORRELATION(QSfHP, YHP)
corry(15)= CORRELATION(TargetHP, YHP)

corrfi = CORRELATION(FSHP, NIIHP)
corrysi = CORRELATION(FS, NII)

!!!!!!!!!!!!!!!!!!!!!!!
! Display the Results !
!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,'(1X,A, A,A,I4)', ADVANCE = 'YES') ' ', dynamicfile, ' iteration: ', iteration	
WRITE(*,*) ' '
WRITE(*,*)
WRITE(*,'(1X,4(A,F7.4),A)', ADVANCE = 'NO') ' [KMlowest = ', lowestKM, ', KMhighest = ', highestKM, &
											']    [SMlowest = ', lowestSM, ', SMhighest = ', highestSM,']'
WRITE(*,*)
WRITE(*,*)
WRITE(*,*)

WRITE(*,FMT = '(1X,a)',ADVANCE = 'YES') '    FS     GDP    C      I      TH     K    adjust   NII    Q      S1     M      IS     X     QSf     S*'
WRITE(*,*) 'mean values of:'
WRITE(*,FMT = '(1X,17(F7.3))', ADVANCE = 'YES') avg(1:15)

WRITE(*,*) ' '	
WRITE(*,*) 'percentage standard deviations of:'
WRITE(*,FMT = '(1X, 17(F7.3))', ADVANCE = 'YES') std(1:15)

WRITE(*,*) ' '	
WRITE(*,*) 'relative standard deviations of:'
WRITE(*,FMT = '(1X, 17(F7.3))', ADVANCE = 'YES') std(1:15)/std(2)

WRITE(*,*) ' '	
WRITE(*,*) 'contemporaneous correlations with Output:'
WRITE(*,FMT = '(1X, 17(F7.3))', ADVANCE = 'YES') corry(1:15)
WRITE(*,*) ' '
WRITE(*,'(1X,2(A,F7.4))', ADVANCE = 'YES') 'filtered corr(FS,NII) ', corrfi, '   unfiltered correlation: ', corrysi
WRITE(*,*) ' ' 

WRITE(*,*) ' ' 
WRITE(*,*) ' ' 

WRITE(*,*) ' mean, std., maximum and minimum adjustment rates and correlation with GDPsim '
WRITE(*,*) 
WRITE(*,'(1X,A, 5F10.4)') ' using fraction of firms adjusting:  ', adjustmoments(1,1:5)
WRITE(*,'(1X,A, 5F10.4)') ' using fraction of imbalance closed: ', adjustmoments(2,1:5)
WRITE(*,*) ' ' 
WRITE(*,*) ' ' 

WRITE(*,*) 'well finished then, press a key ' 
READ(*,*) 

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!
!						!	
!   Standard Deviation	!	
!						!	
!!!!!!!!!!!!!!!!!!!!!!!!!

PURE FUNCTION STDDEV(vector)

! Compute Standard Deviation of a vector, 
IMPLICIT NONE
INTEGER(ik):: n
REAL(rk):: vector(:), mean, STDDEV, num
REAL(rk), ALLOCATABLE:: vectoruse(:)
INTENT(IN):: vector

n = SIZE(vector); ALLOCATE(vectoruse(n))

num = DBLE(n)

mean = SUM(vector)/num
vectoruse = (vector - mean)**2.0_rk

! num = DBLE(n - 1_ik) ! Hogg and Craig (1978), 4th edition, page 124.

STDDEV = SUM(vectoruse)/num
STDDEV = DSQRT(STDDEV)

END FUNCTION STDDEV

!!!!!!!!!!!!!!!!!
!				!
!	COVARIANCE	!
!				!
!!!!!!!!!!!!!!!!!

FUNCTION CORRELATION(vec1, vec2, switch)

! Covariance or correlation between two vectors.  Note this program must recalculate the 
! standard deviations of vec1 and vec2 without using STDDEV.  Internal subroutines 
! cannot access other internal subroutines.  

IMPLICIT NONE
OPTIONAL:: switch
INTEGER(ik):: n1, n2, switch
REAL(rk):: vec1(:), vec2(:), mean1, mean2, std1, std2, num, CORRELATION !, STTDEV
REAL(rk), ALLOCATABLE:: vec(:)
INTENT(IN):: vec1, vec2, switch

n1 = SIZE(vec1)
n2 = SIZE(vec2)

IF (n1.NE.n2) THEN
	WRITE(*,*) ' You bungler, I cannot calculate a correlation between different length series '	
	READ(*,*)
	RETURN
ELSE 
	ALLOCATE(vec(n1))
END IF

num = DBLE(n1)


mean1 = SUM(vec1)/num
mean2 = SUM(vec2)/num

vec = (vec1 - mean1)*(vec2 - mean2)

CORRELATION = SUM(vec)/num		! This is the covariance

IF (PRESENT(switch)) THEN
	IF (switch.EQ.0_ik) THEN	! return covariance not correlation coefficient
		RETURN
	END IF	
END IF

! Compute standard deviations
vec = (vec1 - mean1)**2.0_rk
std1 = SUM(vec)/num
std1 = DSQRT(std1)

vec = (vec2 - mean2)**2.0_rk
std2 = SUM(vec)/num
std2 = DSQRT(std2)

CORRELATION = CORRELATION/(std1*std2)	! correlation coefficient is the default

END FUNCTION CORRELATION

END PROGRAM OutputInv


SUBROUTINE HPfilterS(N,Y,lambda, tau)	! Sparse version, 01.12.2001
 
USE KINDSET								! All matrices are 5 columns or less,	
										! instead of N.  This allows for sparse storage of 	
IMPLICIT NONE							! large simulations.  Any element M(i,j) is now 
										! referenced as M(i,rowc + [j-i]).  Conceptually, 
INTEGER(ik):: row, rowc					! rowc always equals row.  However, for compact 
INTEGER(ik), INTENT(IN):: N				! storage we set rowc = 0.   This is an effective
REAL(rk), INTENT(IN):: Y(N), lambda		! device as most matrix references did involve 
REAL(rk), INTENT(OUT):: tau(N)			! M(row, row + l) where l = -2,-1,0,1, or 2.
REAL(rk):: G(N, -2:2), L(N,-2:0), U(N,0:2), X(N)

! HPfilter(N, Y, lambda, tau) operates the Hodrick-Prescott Filter on a (N x 1)
! time series, Y.  lambda is the smoothing parameter and tau is the computed trend
!
!
! Ordinarily lambda = 1600 for quarterly data and lambda = 100 for annual
!
! The first order conditions to locate the trend (TAU) for an T length
! vector of data Y are:  Y = G*TAU.

! Set up the matrix G which applies the filter.
! row = 1, rowc = row, this allows the use of a sparse 
! approach for large simulations.

rowc = 0
G = 0.0_rk
G(1,rowc) = 1.0_rk+lambda
G(1,rowc + 1) = -2.0_rk*lambda
G(1,rowc+2) = lambda
G(2,rowc-1) = -2.0_rk*lambda
G(2,rowc) = 1.0_rk + 5.0_rk*lambda
G(2,rowc+1) = -4.0_rk*lambda
G(2,rowc+2) = lambda
G(N-1,rowc-2) = lambda
G(N-1,rowc-1) = -4.0_rk*lambda
G(N-1,rowc) = 1.0_rk + 5.0_rk*lambda
G(N-1,rowc+1) = -2.0_rk*lambda
G(N,rowc-2) = lambda
G(N,rowc-1) = -2.0_rk*lambda
G(N,rowc) = 1.0_rk+lambda

DO row = 3, N - 2
     G(row,rowc-2) = lambda
     G(row,rowc+2) = G(row,rowc-2)
     G(row,rowc-1) = -4.0_rk*lambda
     G(row,rowc+1) = G(row,rowc-1)
     G(row,rowc) = 1.0_rk + 6.0_rk*lambda
END DO

! LU factorisation of quindiagonal matrix G, algorithm in Notes.  
CALL LU5(N, G,L,U)

! Having LU factorized the matrix G, we now solve for the terms tau in 
! G*tau = y as L(U*tau) = y first solving Lx = y for x

! The forward pass in Lx = y, based on the form for L, see notes of 11.10.2001.
X(1) = Y(1)
X(2) = Y(2) - L(2,rowc-1)*X(1)

DO row = 3,N
	X(row) = Y(row) - L(row,rowc-2)*X(row-2) - L(row,rowc-1)*X(row-1)
END DO
	
! Having retrieved x, the backward pass in U*tau = x.
tau(N) = X(N)/U(N,rowc)
tau(N-1) = (X(N-1) - tau(N)*U(N-1,rowc+1))/U(N-1,rowc)

DO row = N-2,1,-1
	tau(row) = (X(row) - U(row,rowc+2)*tau(row+2) - U(row,rowc+1)*tau(row+1))/U(row,rowc)
END DO

CONTAINS

! Algorithm to compute the LU factorisation of a quindiagonal matrix A 
! that has at most 5 non-zero elements along its principal diagonal.  
! Based on the tridiadonal example in Van Loan (2000)

SUBROUTINE LU5(N,G,L,U)

INTEGER(ik):: row, rowc
INTEGER(ik), INTENT(IN):: N
REAL(rk), INTENT(IN):: G(N,-2:2)
REAL(rk), INTENT(OUT):: L(N,-2:0), U(N,0:2)

rowc = 0
! initialize sparse matrices L and U
L = 0.0_rk
U = 0.0_rk

L(1,rowc) = 1.0_rk
U(1,rowc:rowc+2) = (/G(1,rowc), G(1,rowc+1), G(1,rowc+2)/)

L(2,rowc) = 1.0_rk
L(2,rowc-1) = G(2,rowc-1)/U(1,rowc)
L(3,rowc-2) = G(3,rowc-2)/U(1,rowc)

U(2,rowc) = G(2,rowc) - L(2,rowc-1)*U(1,rowc+1)
U(2,rowc+2) = G(2,rowc+2)

DO row = 3, N-2				
	
	! The unit principal diagonal term for L
	L(row,rowc) = 1.0_rk
	U(row,rowc+2) = G(row,rowc+2)

	! The 2 L elements related to row (we're really going along a column for each case)
	L(row, rowc - 1) = (G(row,rowc-1) - L(row,rowc-2)*U(row-2,rowc+1))/U(row-1,rowc)
	L(row+1,rowc-2) = G(row+1,rowc-2)/U(row-1,rowc)

	! The 2 U elements related to row
	U(row-1,rowc+1) = G(row-1,rowc+1) - L(row-1,rowc-1)*U(row-2,rowc+2)
	U(row,rowc) = G(row,rowc) - L(row,rowc-2)*U(row-2,rowc+2) - L(row,rowc-1)*U(row-1,rowc+1)

END DO

! row N-1, exclusion for lack of upper diagonal in U, no U(row, row+2) term
row = N-1
L(row,rowc) = 1.0_rk
L(row,rowc-1) = (G(row,rowc-1) - L(row,rowc-2)*U(row-2,rowc+1))/U(row-1,rowc)
L(row+1,rowc-2) = G(row+1,rowc-2)/U(row-1,rowc)
U(row-1,rowc+1) = G(row-1,rowc+1) - L(row-1,rowc-1)*U(row-2,rowc+2)
U(row,rowc) = G(row,rowc) - L(row,rowc-2)*U(row-2,rowc+2) - L(row,rowc-1)*U(row-1,rowc+1)

! row N, exclusion for lack of L(row+1,row), U(row, row+1) and U(row, row+2) terms
row = N
L(row,rowc) = 1.0_rk
L(row,rowc-1) = (G(row,rowc-1) - L(row,rowc-2)*U(row-2,rowc+1))/U(row-1,rowc)
U(row-1,rowc+1) = G(row-1,rowc+1) - L(row-1,rowc-1)*U(row-2,rowc+2)
U(row,rowc) = G(row,rowc) - L(row,rowc-2)*U(row-2,rowc+2) - L(row,rowc-1)*U(row-1,rowc+1)

END SUBROUTINE LU5

END SUBROUTINE HPFILTERS