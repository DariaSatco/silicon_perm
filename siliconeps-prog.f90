program siliconeps

use dopcasimir
use spline
use silicon
use quadpack

implicit none

integer, parameter:: n=483 !number of rows in the file
real(8):: matrixSil(n,3), epsteorSil(100)
real(8):: x(100) !frequency vector
!variables needed fo Kr-Kr integral in cspint subroutine
real(8)::funval(n),Spcoef(3,n), exint(n), work(n)
!variables for Kr-K integral in qags subroutine
real(4):: abserr,abserr2
integer(4)::neval,ier,neval2,ier2
!integral-results in Kr-Kr
real(8):: integral1,integral2,integral3
!Silicon permittivity
real(8), allocatable:: epsSil(:)
!parameters for article formula
real(8):: eps0, epsinf, w0
!counters
integer:: i,j,k
real(8):: splineint

open(unit=10, file='resultSi.txt', status='old')
!read matrix from file
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

200 format(' ', e15.3, f15.7, f15.7)

do i=1,n
read(10,*), (matrixSil(n+1-i,j),j=1,3)
end do

!do i=1,n
!write(*,200), (matrixSil(i,j),j=1,3)
!end do


!Kramers-Kronig for 'freq'-vector

allocate(epsSil(100))
freq=20_8

	do i=1,100
	
	!integrates the oscillator-function
	call qags(KramersKr, 0.0, sngl(matrixSil(1,1)), 1.0e-5, 1.0e-3, sngl(integral1), abserr, neval,ier);
	
	!builds vector to integrate from the Palik data
		do k=1,n
		funval(k)=matrixSil(k,3)*matrixSil(k,1)/(matrixSil(k,1)**2+freq**2)
		end do
	!integrates the vector of data
 	call cspint(n, matrixSil(1:483,1), funval, matrixSil(1,1), matrixSil(n,1), Spcoef, exint, work, integral2)

	!integrates the oscillator-function
	call qags(KramersKr, sngl(matrixSil(n,1)), sngl(matrixSil(n,1)*1e4_8), 1.0e-5, 1.0e-3, sngl(integral3), abserr2, neval2, ier2);

	
	!Kr-Kr formula for permittivity 
	epsSil(i)=(integral1+integral2+integral3)*2/pi+1

	x(i)=freq
	freq=freq*10**(0.2_8)

	end do

!print *, epsSil

!theory - Drude-Loretz model
!parameters from the article
eps0=11.87_8
epsinf=1.035_8
w0=6.6e+15_8

!article formula

do i=1,100
    epsteorSil(i)=epsinf+(eps0-epsinf)*w0**2/(w0**2+x(i)**2);
end do

 call spline_b_val(100, x, epsSil, x(30)+1.0e4_8, splineint)
print *, 'splineint=', splineint, epsSil(30)

open(unit=17, file='silicondata.txt', status='replace')

!write data in file
do i=1,100
write(17,*) x(i), epsSil(i), epsteorSil(i)
end do


close(10)
close(17)
deallocate(epsSil)

end program
