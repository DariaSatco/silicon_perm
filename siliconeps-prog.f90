program siliconeps

use silicon

implicit none

real(8), parameter:: pi=3.14159265359

integer, parameter:: n=483 !number of rows in the file
real(8):: matrixSil(n,3), epsteorSil(100), epsSilLor(100)
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

!article formula

do i=1,100
    epsteorSil(i)=epsSil_art(x(i))
    epsSilLor(i)=oscillators_full(x(i))
end do

open(unit=17, file='silicondata.txt', status='replace')

100 format(e10.2, 3(f15.3))

!write data in file
do i=1,100
write(17,100) x(i), epsSil(i), epsteorSil(i), epsSilLor(i)
end do


close(10)
close(17)
deallocate(epsSil)

end program
