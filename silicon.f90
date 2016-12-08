module silicon

real(8), public:: freq
!approximation of Im (eps)
real(8), dimension(3), parameter::g1= (/39.1, 121.7, 35.6/) !(eV)^2
real(8), dimension(3), parameter::gammam1=(/1.05_8, 0.69_8, 0.35_8/)    !eV
real(8), dimension(3), parameter::wm1=(/5.48_8, 4.19_8, 3.51_8/)   !eV

!approximation of Re (eps)
real(8), dimension(3), parameter::g2= (/69.5, 121.7, 36.9/)  !(eV)^2
real(8), dimension(3), parameter::gammam2=(/1.05_8, 0.69_8, 0.35_8/)    !eV
real(8), dimension(3), parameter::wm2= (/5.48_8, 4.19_8, 3.51_8/)   !eV

contains

!--------------------------------------------------------------------------

function oscillators(x)
!imaginary part of eps(w) in the Lorentz model

real(4):: x, oscillators
real(8), dimension(3):: g, gammam, wm

g=(g1+g2)/2
gammam=(gammam1+gammam2)/2
wm=(wm1+wm2)/2

oscillators=0.0
	do i=1,size(g)
	oscillators=oscillators+g(i)*gammam(i)*x/((wm(i)**2-x**2)**2+gammam(i)**2*x**2)
	end do

end function

!-------------------------------------------------------------------------

function oscillators_full(x)
!eps(iw) in the Lorentz model

real(8):: x, oscillators_full
real(8), dimension(3):: g, gammam, wm

g=(g1+g2)/2
gammam=(gammam1+gammam2)/2
wm=(wm1+wm2)/2

oscillators_full=1.0_8
    do i=1,size(g)
    oscillators_full=oscillators_full+g(i)/(wm(i)**2+x**2+gammam(i)*x)
    end do

end function


!-----------------------------------------------------------------------

function KramersKr(x)

real(4):: x, KramersKr

KramersKr=x/(x**2+freq**2)*oscillators(x)

end function


!---------------------------------------------------------------------

function epsSil_art(x)

real(8):: epsSil_art
real(8):: x
real(8):: eps0, epsinf, w0

!silicon parameters
!theory - Drude-Loretz model
!parameters from the article
eps0=11.87_8
epsinf=1.035_8
w0=4.34_8

epsSil_art=epsinf+(eps0-epsinf)*w0**2/(w0**2+x**2)

end function

!---------------------------------------------------------------------

end module
