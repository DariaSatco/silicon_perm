program siliconeps

use silicon

implicit none

real(8), parameter:: pi=3.14159265359
real(8), parameter :: kb=8.617e-5  ! eV/K kb=1.38e-23 J/K  Boltzman constant
real(8), parameter :: h=6.585e-16  ! eV*s h=1.054e-34J*s  Plank constant

integer, parameter:: n=483 !number of rows in the file
integer, parameter:: nw=100 !number of Matsubara frequencies

real(8):: matrixSil(n,3), epsteorSil(nw), epsSilLor(nw), epsSil(nw)
real(8):: x(nw) !frequency vector
!variables needed fo Kr-Kr integral in cspint subroutine
real(8)::funval(n),Spcoef(3,n), exint(n), work(n)
!integral-results in Kr-Kr
real(8):: integral
real(8):: eV, T
!counters
integer:: i,j,k

print*, 'Process...'
open(unit=10, file='resultSi.txt', status='old')
!read matrix from file
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

open(unit=11, file='resultSi_eV.txt', status='replace')

eV=1.519e15 !rad/s
T=5._8

200 format(3(f15.3))

do i=1,n
read(10,*), (matrixSil(n+1-i,j),j=1,3)
end do

do i=1,n
matrixSil(i,1)=matrixSil(i,1)/eV
write(11,200), (matrixSil(i,j),j=1,3)
end do


!Kramers-Kronig for 'freq'-vector

	do i=1,nw
	
	x(i)=2*pi*kb*T*i/h
    x(i)=x(i)/eV
    freq=x(i)

	!build vector to integrate from the Palik data
		do k=1,n
		funval(k)=matrixSil(k,3)*matrixSil(k,1)/(matrixSil(k,1)**2+freq**2)
		end do
	!integrate the vector of data
 	call cspint(n, matrixSil(1:483,1), funval, matrixSil(1,1), matrixSil(n,1), Spcoef, exint, work, integral)
	
	!Kr-Kr formula for permittivity 
	epsSil(i)=integral*2/pi+1

	end do

do i=1,nw
    epsteorSil(i)=epsSil_art(x(i))
    epsSilLor(i)=oscillators_full(x(i))
end do

open(unit=17, file='silicondata.txt', status='replace')

100 format(4(f15.3))

!write data in file
do i=1,nw
write(17,100) x(i), epsSil(i), epsteorSil(i), epsSilLor(i)
end do


close(10)
close(11)
close(17)

print*, 'Done!'

end program
