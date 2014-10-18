Program LinearEquations
! solving the matrix equation A*x=b using LAPACK
!!!!!!! A must be NXN !!!!!!!

Implicit none

! declarations, notice single precision
Real :: A(3,3), b(3), btemp(3), x(3), AT(3,3)
integer i, pivot(3), ok

! define matrix A

! first row
A(1,1) = 1.
A(1,2) = 2.  
A(1,3) = 3.

! second row
A(2,1) = 5.
A(2,2) = 9.  
A(2,3) = 0.

! third row
A(3,1) = 1.
A(3,2) = 4.  
A(3,3) = 732.


AT = TRANSPOSE(A)

! define vector b, make b a matrix and you can solve multiple
! equations with the same A but different b
b(1) = 1.
b(2) = 2.
b(3) = 435.

btemp(1) = 1.
btemp(2) = 2.
btemp(3) = 435.

! find the solution using the LAPACK routine SGESV
call GESV(3, 1, A, 3, pivot, b, 3, ok)



! parameters in the order as they appear in the function call
!    order of matrix A, number of right hand sides (b), matrix A,
!    leading dimension of A, array that records pivoting, 
!    result vector b on entry, x on exit, leading dimension of b
!    return value 
	
!multiply
x = MATMUL(b,AT)
	
! print the vector x
do i = 1, 3
  write(*,*) b(i)
end do
end


 ! run with the following commands:
 ! gfortran -o lapak_ex lapak_ex.f90 -L/usr/lib -llapack -lblas
 ! ./lapak_ex
