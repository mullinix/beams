PROGRAM MAT_TEST

Implicit none
EXTERNAL DGESV
! declarations, notice single precision
INTEGER N
PARAMETER (N=256)
INTEGER I, pivot(N), ok
INTEGER T1,T2, CLOCK_RATE, CLOCK_MAX
REAL TIMES

DOUBLE PRECISION A(N,N)
DOUBLE PRECISION B(N,1)
DOUBLE PRECISION C(N,1)
DOUBLE PRECISION X(N,1)

CALL RANDOM_NUMBER(A)
CALL RANDOM_NUMBER(B)
C=B
CALL SYSTEM_CLOCK ( T1 , CLOCK_RATE, CLOCK_MAX ) 
DO I=1,1000
! find the solution using the LAPACK routine SGESV
CALL DGESV(N, 1, A, N, pivot, B, N, ok)
X=MATMUL(A,B)
B=C
END DO

CALL SYSTEM_CLOCK ( T2 , CLOCK_RATE, CLOCK_MAX )
TIMES = REAL(T2 - T1)/REAL(CLOCK_RATE)
PRINT*, TIMES
! parameters in the order as they appear in the function call
!    order of matrix A, 
!	 number of right hand sides (B), 
!	 matrix A,
!    leading dimension of A, 
!	 array that records pivoting, 
!    result vector b on entry, 
!	 x on exit, 
!    leading dimension of b
!    return value 

END PROGRAM
