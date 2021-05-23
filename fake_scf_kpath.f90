!programa para gerar o kpath para fake_scf de funcionais hibridos

!gfortran fake_scf_kpath.f90 -o f_scf_kpath.x
! f_scf_kpath.x <input_fake_scf.dat> output_fake_scf.dat

program main

	implicit none

	integer :: i,j
	integer :: npt
	integer :: ntr
	double precision,allocatable,dimension(:,:) :: kpt
	double precision :: x,y,z

	read(*,*) npt
	read(*,*) ntr

	allocate(kpt(npt,3))

	do i=1,npt
	 read(*,*) kpt(i,1),kpt(i,2),kpt(i,3)
	end do	

	do i=1,npt-1

		do j=1,ntr

		x=kpt(i,1)+(kpt(i+1,1)-kpt(i,1))*(dble(j-1)/dble(ntr))
		y=kpt(i,2)+(kpt(i+1,2)-kpt(i,2))*(dble(j-1)/dble(ntr))
		z=kpt(i,3)+(kpt(i+1,3)-kpt(i,3))*(dble(j-1)/dble(ntr))

		write(*,*) x,y,z,0
	
		end do

	        write(*,*) kpt(i+1,1),kpt(i+1,2),kpt(i+1,3),0

	end do


 


end program main
