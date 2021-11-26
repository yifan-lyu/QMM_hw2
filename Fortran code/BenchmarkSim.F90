SUBROUTINE BenchmarkSim(Nz, Z, PI, rknotsk, alpha, thetam, thetan, &
                        Atech, delta, beta, precisiond, Apref, knotsk, cskd, kbound, tauchenscale1, tauchenscale2, &
                        datafile, shockfile, simlength, again, ZIsim, Qsim, Psim, Ksim)

! Simulate the benchmark model of Khan and Thomas (2007).

USE KINDSET

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

end interface



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                     !
!   Subroutine BenchmarkSim.f90 Type Declarations     !
!                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


INTEGER(ik):: Nz, rknotsk, rknotskpoly, rknotskint, simlength, time, iz, ZIsim(simlength), varnum, &
			  lags, j, again, seedset(2)
REAL(rk):: Z(2, Nz), PI(Nz, Nz), ll, lh, cskd(4,rknotsk-1, Nz), switch, simnum, lambdasmooth, &
		   alpha, thetam, thetan, Atech, delta, beta, precisiond, Apref, &
		   kbound(2), avg(15), std(15), corrY(15), zval, gzval, kval, knotsk(rknotsk), errorval, begin, &
		   finish, realrate, Zsim(2,simlength), Ysim(simlength), Csim(simlength), Isim(simlength), &
		   Ksim(simlength+1), THsim(simlength), Lambdasim(simlength), NMsim(simlength), mur, &
		   NKsim(simlength), Msim(simlength), YHP(simlength), CHP(simlength), IHP(simlength), &
		   THHP(simlength), MHP(simlength), &
		   Z1HP(simlength), Z2HP(simlength), KHP(simlength), Qsim(simlength), GDPsim(simlength), QHP(simlength), &
		   Psim(simlength), QMHP(simlength), MUmultsim(simlength), pmean, &
		   qcheck, mshare, zstats(7), Nshare(simlength), KoverN(simlength), MoverN(simlength), &
		   tauchenscale1, tauchenscale2

REAL(rk), ALLOCATABLE:: TREND(:)

CHARACTER(30):: datafile, shockfile, resultfile

INTENT(IN):: Nz, Z, PI, rknotsk, alpha, thetam, thetan, &
			 Atech, delta, beta, precisiond, Apref, knotsk, cskd, kbound, &
			 tauchenscale1, tauchenscale2, again, simlength, shockfile
INTENT(OUT):: ZIsim, Qsim, Psim, Ksim

EXTERNAL HPFILTERS


!!!!!!!!!!!!!
!			!
!	Setup	!
!			!
!!!!!!!!!!!!!

rknotskint = rknotsk - 2_ik; rknotskpoly = rknotsk - 1_ik

CALL CPU_TIME(begin)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!								!
!	Simulate discrete shock		!
!								!		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (again.eq.-1) THEN	! resolve a simulation of total factor productivity

	DO 

		WRITE(*,FMT = '(1X,A)', ADVANCE = 'NO') ' enter existing shock file: '
		READ(*,*) datafile	
	
		OPEN (UNIT = 57, FILE = datafile, STATUS = 'OLD', ACTION = 'READ', IOSTAT = lags, ERR = 100)

		! The subroutine goes to line 100, right below, if there's a file opening error.  
		100	IF (lags.NE.0_ik) THEN 
				WRITE(*,*) ' Hey, I think I could not locate ', datafile
				WRITE(*,*) ' I got a file I/O error ', lags
			ELSE
				EXIT
			END IF

	END DO

	READ(57,*) j

	DO time = 1,simlength
		READ(57,*) zIsim(time), Qsim(time), Psim(time)
	END DO

	DO time = 1, simlength + 1
		READ(57,*) Ksim(time)
	END DO

	CLOSE(57)
    
    DO time = 1, simlength
        Zsim(1,time) = Z(1,ZIsim(time))
        Zsim(2,time) = Z(2,ZIsim(time))
    END DO

ELSE

	CALL ZSHOCK2(Nz, Z, PI, simlength, seedset, Zsim, ZIsim)

END IF 

zstats(1) = SUM(DLOG(Zsim(1,1:simlength)))/DBLE(simlength)
zstats(2) = STDDEV(DLOG(Zsim(1,1:simlength)))
zstats(3) = CORRELATION(DLOG(Zsim(1,2:simlength)), DLOG(Zsim(1,1:simlength-1)))
zstats(4) = SUM(DLOG(Zsim(2,1:simlength)))/DBLE(simlength)
zstats(5) = STDDEV(DLOG(Zsim(2,1:simlength)))
zstats(6) = CORRELATION(DLOG(Zsim(2,2:simlength)), DLOG(Zsim(2,1:simlength-1)))
zstats(7) = CORRELATION(DLOG(Zsim(1,1:simlength)), DLOG(Zsim(2,1:simlength)))

WRITE(*,'(1X,3(A,F8.4))', ADVANCE = 'YES') ' total factor productivity log z1 mean: ', &
		 zstats(1), ' std. : ', zstats(2), ' autocorr. : ', zstats(3)
WRITE(*,'(1X,3(A,F8.4))', ADVANCE = 'YES') ' total factor productivity log z2 mean: ', &
		 zstats(4), ' std. : ', zstats(5), ' autocorr. : ', zstats(6)		 
WRITE(*,'(1X,A,F8.4)', ADVANCE = 'YES')    ' correlation between log z1 and log z2: ', zstats(7)		 


DO

	IF (again.EQ.-1) THEN	! using an existing simulation with a pre-existing initial condition
		CONTINUE
	ELSE					! this is a new simulation and I need an initial condition
		WRITE(*,FMT = '(1X,A)', ADVANCE = 'NO') ' enter initial state (k): '
		READ(*,*)  Ksim(1)	
	END IF

	IF (Ksim(1).LT.kbound(1).OR.Ksim(1).GT.kbound(2)) THEN
		WRITE(*,*) ' You violated the knots! '
		WRITE(*,*) kbound
	ELSE
		EXIT
	END IF

END DO


! Standard file somedata will deliver the actual shock values and some parameters to OutputInv file
OPEN(UNIT = 35, FILE = 'somedata', ACTION = 'READWRITE', STATUS = 'REPLACE')

DO time = 1,simlength
	WRITE(35,*) (/Zsim(1,time), Zsim(2,time)/)
END DO
WRITE(35,*) Apref, beta

CLOSE(35)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!							!
!	Simulate the economy	!
!							!		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO time = 1, simlength
	
	iz = ZIsim(time)
	zval = Z(1,iz)		! exogenous state determined by simulation program
	gzval = Z(2,iz)     ! second exogenous state determined by simulation program
	kval = Ksim(time)	! predetermined state
	
	CALL lambdabound(zval, gzval, kval, kbound(1), alpha, Apref, delta, thetam, thetan, Atech, ll, precisiond)
	CALL lambdabound(zval, gzval, kval, kbound(2), alpha, Apref, delta, thetam, thetan, Atech, lh, precisiond)

	CALL DR(alpha, thetam, thetan, Atech, delta, Apref, rknotsk, knotsk, cskd(:,:,iz), zval, gzval, &
	        kval, Lambdasim(time), MUmultsim(time), Csim(time), Msim(time), NMsim(time), NKsim(time), &
	        Ksim(time+1), Ysim(time), switch, errorval, precisiond, ll, lh)
			  		
	THsim(time) = NMsim(time) + NKsim(time)
	Qsim(time) = thetam*Ysim(time)/Msim(time)
	Isim(time) = Ksim(time+1) - (1.0_rk - delta)*Ksim(time)
	GDPsim(time) = Ysim(time)
	Psim(time) = 1.0_rk/Csim(time)

	WRITE(*,FMT = '(1X,A,I5,A, 3(F10.4), A, 6(F8.3))', ADVANCE = 'YES') &
	' t = ', time, ' (z,gz,k) = ',zval, gzval, kval, '  |  y, c,i,m, q, th = ', &
	GDPsim(time), Csim(time), Isim(time), Msim(time), Qsim(time), THsim(time)

END DO


simnum = DBLE(simlength)
avg(1) = SUM(Zsim(1,1:simlength))/simnum
avg(2) = SUM(GDPsim)/simnum
avg(3) = SUM(Csim)/simnum
avg(4) = SUM(Isim)/simnum
avg(5) = SUM(THsim)/simnum
avg(6) = SUM(Ksim)/(simnum+1.0_rk)
avg(7) = SUM(Msim)/simnum
avg(8) = SUM(Qsim)/simnum
avg(9) = SUM(Qsim*Msim)/simnum
avg(10) = SUM(Zsim(2,1:simlength))/simnum

realrate = SUM(Csim(2:simlength)/(beta*Csim(1:simlength)) - 1.0_rk)/(simnum - 1.0_rk)
mshare = SUM((Qsim*Msim)/GDPsim)/simnum

! realrate = ((1.0_rk + realrate)**4.0_rk - 1.0_rk)*100.0_rk
pmean = SUM(1.0_rk/Csim)/simnum

! Certainty value of q if p were constant at pmean and z were at avg(1).
qcheck = ((1.0_rk - beta*(1.0_rk - delta))/(beta*alpha))**alpha
qcheck = qcheck*((Apref/(1.0_rk - alpha))**(1.0_rk - alpha))
qcheck = (qcheck)*(pmean**(alpha - 1.0_rk))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute Hodrick-Prescott Filtered Statistics !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

lambdasmooth = 1600.0_rk

ALLOCATE(TREND(simlength))

YHP = DLOG(GDPsim)
CALL HPFILTERS(simlength, YHP, lambdasmooth, TREND)
YHP = YHP - TREND
CHP = DLOG(Csim)
CALL HPFILTERS(simlength, CHP, lambdasmooth, TREND)
CHP = CHP - TREND
IHP = Isim					! separate method of detrending investment in case I(t)<0
CALL HPFILTERS(simlength, IHP, lambdasmooth, TREND)
IHP = (IHP - TREND)/TREND

Z1HP = DLOG(Zsim(1,1:simlength))
CALL HPFILTERS(simlength, Z1HP, lambdasmooth, TREND)
Z1HP = Z1HP - TREND
Z2HP = DLOG(Zsim(2,1:simlength))
CALL HPFILTERS(simlength, Z2HP, lambdasmooth, TREND)
Z2HP = Z2HP - TREND
THHP = DLOG(THsim)
CALL HPFILTERS(simlength, THHP, lambdasmooth, TREND)
THHP = THHP - TREND
MHP = DLOG(Msim)
CALL HPFILTERS(simlength, MHP, lambdasmooth, TREND)
MHP = MHP - TREND
QHP = DLOG(Qsim)
CALL HPFILTERS(simlength, QHP, lambdasmooth, TREND)
QHP = QHP - TREND
QMHP = DLOG(Qsim*Msim)
CALL HPFILTERS(simlength, QMHP, lambdasmooth, TREND)
QMHP = QMHP - TREND

DEALLOCATE(TREND)

ALLOCATE(TREND(simlength+1))

CALL HPFILTERS(simlength+1, DLOG(Ksim), lambdasmooth, TREND)
KHP = DLOG(Ksim) - TREND


DEALLOCATE(TREND)

std(1) = STDDEV(Z1HP)		!Standard Deviations of each Series
std(2) = STDDEV(YHP)
std(3) = STDDEV(CHP)
std(4) = STDDEV(IHP)
std(5) = STDDEV(THHP)
std(6) = STDDEV(KHP)
std(7) = STDDEV(MHP)
std(8) = STDDEV(QHP)
std(9) = STDDEV(QMHP)
std(10) = STDDEV(Z2HP)
std = std*100.0_rk	! report percentage standard deviations as percentages


corry(1)= CORRELATION(Z1HP, YHP)		  !Cross Correlations with Output
corry(2)= CORRELATION(YHP, YHP)
corry(3)= CORRELATION(CHP, YHP)
corry(4)= CORRELATION(IHP, YHP)
corry(5)= CORRELATION(THHP, YHP)
corry(6)= CORRELATION(KHP(1_ik:simlength), YHP)
corry(7)= CORRELATION(MHP, YHP)
corry(8)= CORRELATION(QHP, YHP)
corry(9)= CORRELATION(QMHP, YHP)
corry(10)= CORRELATION(Z2HP, YHP)


!!!!!!!!!!!!!!!!!!!!!!!
! Display the Results !
!!!!!!!!!!!!!!!!!!!!!!!

CALL CPU_TIME(finish); finish = finish - begin

WRITE(*,*) ' ' 
WRITE(*,'(1X,3(A,F5.3))', ADVANCE = 'NO') ' alpha = ', alpha, ' thetam = ', thetam, ' thetan = ', thetan
WRITE(*,'(1x,a,2f10.6)' ) ' tauchenscales = ', tauchenscale1, tauchenscale2
WRITE(*,'(1X, A,F8.4,A,F8.4, A, A, A)', ADVANCE = 'YES') ' K in (', MINVAL(Ksim), ', ', MAXVAL(Ksim), ') ', &
' input: ', datafile
WRITE(*,*) ' ' 

WRITE(*,'(1X,3(A,F8.4), A, 2I8)', ADVANCE = 'YES') ' log z1 mean: ', &
		 zstats(1), ' std. : ', zstats(2), ' autocorr. : ', zstats(3), '    seed = ', seedset 
WRITE(*,'(1X,3(A,F8.4))', ADVANCE = 'YES') ' log z2 mean: ', &
		 zstats(4), ' std. : ', zstats(5), ' autocorr. : ', zstats(6)
WRITE(*,'(1X,A,F8.4)', ADVANCE = 'YES')    ' correlation between log z1 and log z2: ', zstats(7)			 
WRITE(*,*) ' '

WRITE(*,*) ' ' 
WRITE(*,FMT = '(1X,a)',ADVANCE = 'YES') '     Z1     Z2     GDP      C       I      TH      K     M       Q       QM '
WRITE(*,*) 'Mean Values of:'
WRITE(*,FMT = '(1X,15(F8.3))', ADVANCE = 'YES') (/avg(1), avg(10), avg(2:9)/)

WRITE(*,*) ' '	
WRITE(*,*) 'percentage standard deviations of:'
WRITE(*,FMT = '(1X,15(F8.3))', ADVANCE = 'YES') (/std(1), std(10), std(2:9)/)

WRITE(*,*) ' '	
WRITE(*,*) 'relative standard deviations of:'
WRITE(*,FMT = '(1X,15(F8.3))', ADVANCE = 'YES') (/std(1)/std(2), std(10)/std(2), std(2:9)/std(2)/)

WRITE(*,*) ' '	
WRITE(*,*) 'contemporaneous correlations with Output:'
WRITE(*,FMT = '(1X,15(F8.3))', ADVANCE = 'YES') (/corry(1), corry(10), corry(2:9)/)
WRITE(*,*) ' '

WRITE(*,'(1X, 5(A, F6.3))', ADVANCE = 'YES') &
' avg. K/Y = ', avg(8)/avg(2), ' (annual= ', avg(8)/(avg(2)*4.0_rk),') quarterly rate = ', realrate, ' mean p = ', pmean, ' certainty q = ', qcheck
WRITE(*,*) ' ' 

WRITE(*,*) ' ' 

resultfile = shockfile; j = LEN_TRIM(datafile)

WRITE(*,'(1X,F6.3, A)', ADVANCE = 'NO')  finish, ' second simulation, press enter '


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!							!!
!!							!!
!!	 Internal Subroutines	!!
!!							!!
!!							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!											!
!	Simulate Discrete Stochastic Process	!
!											!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
SUBROUTINE ZSHOCK2(Nz, Z, PI, T, seedset, Zsim, ZIsim)

INTEGER(ik):: Nz, T, time, col, row, ZIsim(T), seedsize, seedset(2)
INTEGER(ik), ALLOCATABLE:: minplace(:), SEED(:)
REAL(rk):: Z(2,Nz), PI(Nz,Nz), Zsim(2,T), RV(T-1), PISUM(Nz, Nz), PID(Nz)

INTENT(IN):: Nz, Z, PI, T
INTENT(OUT):: Zsim, ZIsim, seedset

! automatic integer truncation
ZIsim(1) = (Nz + 1_ik)/2_ik
Zsim(1,1) = Z(1,ZIsim(1))
Zsim(2,1) = Z(2,ZIsim(1))

PISUM = 0.0_rk
PISUM(:,1) = PI(:,1)

DO col = 2, Nz

	PISUM(:,col) = PISUM(:,col - 1) + PI(:,col)

END DO

CALL RANDOM_SEED (SIZE = seedsize)  ! size of seed
ALLOCATE(SEED(seedsize))
SEED = (/345667,39987/)			
seedset = SEED

CALL RANDOM_SEED (PUT = SEED)   

CALL RANDOM_NUMBER(RV)
ALLOCATE(minplace(1))

DO time = 2, T

	! the current Z values index position is row
	row = ZIsim(time - 1)
	! subtract scalar from vector implies subtracting it from each element
	PID = PISUM(row,:) - RV(time-1)
	
	! find least upper bound in PISUM for RV
	minplace = MINLOC(PID, MASK = PID.GT.0.0_rk)
	ZIsim(time) = minplace(1)
	Zsim(1,time) = Z(1,minplace(1))
	Zsim(2,time) = Z(2,minplace(1))
END DO

DEALLOCATE(minplace)

END SUBROUTINE ZSHOCK2


!!!!!!!!!!!!!!!!!!!!!!!!!!
!						!	
!   Standard Deviation	!	
!						!	
!!!!!!!!!!!!!!!!!!!!!!!!!

PURE FUNCTION STDDEV(vector)

! Compute Standard Deviation of a vector 
IMPLICIT NONE
INTEGER(ik):: n
REAL(rk):: vector(:), mean, STDDEV, num
REAL(rk), ALLOCATABLE:: vectoruse(:)

INTENT(IN):: vector

n = SIZE(vector); ALLOCATE(vectoruse(n))

num = DBLE(n)

mean = SUM(vector)/num
vectoruse = (vector - mean)**2.0_rk

STDDEV = SUM(vectoruse)/num
STDDEV = DSQRT(STDDEV)

END FUNCTION STDDEV

!!!!!!!!!!!!!!!!!
!				!
!	COVARIANCE	!
!				!
!!!!!!!!!!!!!!!!!

FUNCTION CORRELATION(vec1, vec2, switch)

! Covariance or correlation between two vectors

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

END SUBROUTINE BenchmarkSim


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