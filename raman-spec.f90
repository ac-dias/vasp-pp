program main

	implicit none

	integer :: i,j	
	integer :: natoms 
	integer :: nmodes
	double precision,parameter :: pi=acos(-1.)

	double precision :: freqmin
	double precision :: freqmax
	integer,parameter :: nfreq= 50001
	double precision :: sme
	

	double precision :: gaussian	
	double precision :: freq,soma	
	double precision,allocatable,dimension(:,:) :: ramandata
	double precision :: aux
	character(len=2) :: caux
	
	!variaveis getarg
	
	character*30 arg
  	integer ios
	
      CALL getarg(1, arg)
		read(arg,*,iostat=ios) natoms

      CALL getarg(2, arg)
		read(arg,*,iostat=ios) freqmin
		
      CALL getarg(3, arg)
		read(arg,*,iostat=ios) freqmax
		
      CALL getarg(4, arg)
		read(arg,*,iostat=ios) sme					


	nmodes = 3*natoms
	
	allocate(ramandata(nmodes-3,2))
	
	read(*,*) caux
	
	do i=1,nmodes-3
	
		read(*,*) aux,ramandata(i,1),aux,aux,ramandata(i,2)
	
	end do
	
	
	do i=1,nfreq
	
		freq= freqmin + (freqmax-freqmin)*((i-1.0)/(nfreq-1.0))
		
		soma = 0.0
		
		do j=1,nmodes
		
			soma=soma+gaussian(freq-ramandata(j,1),sme)*ramandata(j,2)
		
		end do

		write(*,*) freq,(soma)
	
	end do	

end program main


function gaussian(deltaen,sme)

	implicit none
	
	double precision :: gaussian
	double precision :: deltaen,sme
	double precision,parameter :: pi=acos(-1.)
	
	double precision :: norm,deno,aux1	
	
	norm = 1.0D0/(sme*DSQRT(2.0D0*pi))
	deno = 2.0D0*sme*sme
	aux1 = -1.0D0*((deltaen)*(deltaen))/deno
	
	gaussian = norm*dexp(aux1)	

end function gaussian
