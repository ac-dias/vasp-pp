program main

	implicit none

	integer :: i,j	
	integer :: natoms,erro 
	integer :: nmodes
	double precision,parameter :: pi=acos(-1.)

	double precision :: freqmin
	double precision :: freqmax
	integer,parameter :: nfreq= 50001
	double precision :: sme
	

	double precision :: gaussian	
	double precision :: freq,soma	
	double precision,allocatable,dimension(:,:) :: specdata
	double precision :: aux
	character(len=2) :: caux,tipo
	
		
	
	!variaveis getarg
	
	
	write(*,*) "Select 'R' for Raman or 'IR' for Infrared Spectrum"
	read(*,*) tipo
	
	
	if ( (tipo .eq. 'R' ) .or. (tipo .eq. 'IR')  ) then
	
	
	continue
	
	else
	
		write(*,*) "This option is not available"
		STOP
	
	end if
	
	write(*,*) "Number of atoms in the unit cell"
	read(*,*) natoms
	
	write(*,*) "Minimum frequency in cm^-1"
	read(*,*) freqmin
	
	write(*,*) "Maximum frequency in cm^-1"
	read(*,*) freqmax
	
	write(*,*) "Gaussian Smearing in cm^-1"
	read(*,*) sme
	

	
		
	nmodes = 3*natoms
	
	
	allocate(specdata(nmodes-3,2))
	
	select case (tipo)
	
	case('R')
	
	 !INPUT : raman spectrum input
	OPEN(UNIT=100, FILE= "./vasp_raman.dat",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening raman spec input file"
	! OUTPUT : raman spectrum output
	OPEN(UNIT=200, FILE= "raman_spec.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening raman spec output file"
    	
    	read(100,*) caux
	
	do i=1,nmodes-3
	
		read(100,*) aux,specdata(i,1),aux,aux,specdata(i,2)
	
	end do
	
	case('IR')

	! INPUT : ir spectrum input
	OPEN(UNIT=100, FILE= "./intensities/results/exact.res.txt",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening ir spec input file"	
	! OUTPUT : ir spectrum output
	OPEN(UNIT=200, FILE= "ir_spec.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening ir spec output file"
    	
    	do i=1,nmodes-3
	
		read(100,*) aux,specdata(i,1),specdata(i,2)
	
	end do
	
	case default
	 stop
	
	end select
	
	
	

	
	
	do i=1,nfreq
	
		freq= freqmin + (freqmax-freqmin)*((i-1.0)/(nfreq-1.0))
		
		soma = 0.0
		
		do j=1,nmodes
		
			soma=soma+gaussian(freq-specdata(j,1),sme)*specdata(j,2)
		
		end do

		write(200,*) freq,soma
	
	end do	
	
	
	deallocate (specdata)
	close(100)
	close(200)

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
