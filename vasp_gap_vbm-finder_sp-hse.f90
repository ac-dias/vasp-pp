!programa para achar o bandgap
!gfortran vasp_gap_vbm-finder_sp-hse.f90 -o gapf_sp_h.x
!gapf_sp_h.x 30 , 30 numero de pontos k com peso diferente de 0

program main

	implicit none

	integer :: nbands,nkpts
	integer :: erro

	integer :: i,j

	double precision :: flag

	double precision,allocatable,dimension(:,:) :: kpt
	integer,allocatable,dimension(:) :: nocp,nocp2

	double precision,allocatable,dimension(:,:,:) :: energy

	double precision :: gapdir,cbm,vbm
	integer :: cbmkpt,vbmkpt,gapdirkpt
	double precision :: gaux,vaux,caux

	double precision :: gapdir2,cbm2,vbm2
	integer :: cbmkpt2,vbmkpt2,gapdirkpt2
	double precision :: gaux2,vaux2,caux2

	integer ::  gapdirkptf
	double precision :: gapfaux
	

	double precision :: gapdirf,cbmf,vbmf
	double precision :: cbmdir,vbmdir

	double precision,parameter :: ocptol=0.05

	!variaveis getarg

  	character*30 arg
  	integer ios

	! INPUT : lendo os parametros do modelo de tight-binding
	OPEN(UNIT=100, FILE= "EIGENVAL",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de entrada 1"

	! OUTPUT : lendo os parametros do modelo de tight-binding
	OPEN(UNIT=200, FILE= "gap_vasp_sp.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida 1"

      CALL getarg(1, arg)
		read(arg,*,iostat=ios) nkpts
       !CALL getarg(2, arg)
	!	read(arg,*,iostat=ios) nkpts
 


	!nbands=20
	!nkpts=27



	do i=1,5
		read(100,*)
	end do

	read(100,*) flag,flag,nbands

	allocate(kpt(nkpts,3),nocp(nkpts),energy(nkpts,nbands,4))
	allocate(nocp2(nkpts))

	!lendo os valores de energia nos pontos k

	do i=1,nkpts

		read(100,*) kpt(i,1),kpt(i,2),kpt(i,3),flag

		do j=1,nbands

			read(100,*) flag,energy(i,j,1),energy(i,j,2),energy(i,j,3),energy(i,j,4)


		end do

		

	end do

	!termino da leitura

!####################### parte spin 1

	!identificando a primeira banda desocupada de cada ponto k - spin 1

	do i=1,nkpts


		do j=1,nbands

			if (energy(i,j,3) .lt. 1.0) then

				nocp(i)=j


				go to 50

			else 

				continue

			end if

		end do

50	continue		

	end do

	!termino


	gapdir=500
	cbm=500
	vbm=-500

	!calculando os gaps

	!teste metalico

	do i=1,nkpts

		if ( (energy(i,nocp(i),3) .gt. 0.0) .and. (energy(i,nocp(i),3) .lt. 1.0-ocptol)) then

			go to 70

		else 

			continue

		end if

	end do
		
	!fim teste metalico
	cbmkpt=1
	vbmkpt=1
	gapdirkpt=1

	do i=1,nkpts

		if (energy(i,nocp(i),1) .lt. cbm) then

			cbm=energy(i,nocp(i),1)
			cbmkpt=i

		else

			continue
		end if

		if (energy(i,nocp(i)-1,1) .gt. vbm) then

			vbm=energy(i,nocp(i)-1,1)
			vbmkpt=i

		else

			continue
		end if


		gaux=energy(i,nocp(i),1)-energy(i,nocp(i)-1,1)

		if (gaux .lt. gapdir) then

			gapdir=gaux
			gapdirkpt=i
		else

		end if

	end do

	write(200,*) "spin 1"
	write(200,*) "spin 1 gap direto:",real(gapdir)
	write(200,*) "spin 1 gap indireto", real(cbm-vbm)
	write(200,*)
	write(200,*) "spin 1 cbm", real(cbm)
	write(200,*) "spin 1 vbm", real(vbm)
	write(200,*)
	write(200,*) "spin 1 ponto k cbm:",real(kpt(cbmkpt,1)),real(kpt(cbmkpt,2)),real(kpt(cbmkpt,3))
	write(200,*) "spin 1 ponto k vbm:",real(kpt(vbmkpt,1)),real(kpt(vbmkpt,2)),real(kpt(vbmkpt,3))
	write(200,*)
	write(200,*) "spin 1 ponto k gap direto:",real(kpt(gapdirkpt,1)),real(kpt(gapdirkpt,2)),real(kpt(gapdirkpt,3))
	write(200,*) "spin 1 gp dir cb:", real(energy(gapdirkpt,nocp(gapdirkpt),1))
	write(200,*) "spin 1 gp dir vb:",real(energy(gapdirkpt,nocp(gapdirkpt)-1,1))

	go to 60

70	write(200,*) "sistema metalico"


60	continue


!####################### parte spin 2

	!identificando a primeira banda desocupada de cada ponto k - spin 2

	do i=1,nkpts


		do j=1,nbands

			if (energy(i,j,4) .lt. 1.0) then

				nocp2(i)=j


				go to 100

			else 

				continue

			end if

		end do

100	continue		

	end do	

	gapdir2=500
	cbm2=500
	vbm2=-500

	!calculando os gaps

	!teste metalico

	do i=1,nkpts

		if ( (energy(i,nocp2(i),4) .gt. 0.0) .and. (energy(i,nocp2(i),4) .lt. 1.0-ocptol)) then

			go to 140

		else 

			continue

		end if

	end do
		
	!fim teste metalico
	cbmkpt2=1
	vbmkpt2=1
	gapdirkpt2=1

	do i=1,nkpts

		if (energy(i,nocp2(i),2) .lt. cbm2) then

			cbm2=energy(i,nocp2(i),2)
			cbmkpt2=i

		else

			continue
		end if

		if (energy(i,nocp2(i)-1,2) .gt. vbm2) then

			vbm2=energy(i,nocp2(i)-1,2)
			vbmkpt2=i

		else

			continue
		end if


		gaux2=energy(i,nocp2(i),2)-energy(i,nocp2(i)-1,2)

		if (gaux2 .lt. gapdir2) then

			gapdir2=gaux2
			gapdirkpt2=i
		else

		end if

	end do
	
	write(200,*) "                           "
	write(200,*) "spin 2"
	write(200,*) "spin 2 gap direto:",real(gapdir2)
	write(200,*) "spin 2 gap indireto", real(cbm2-vbm2)
	write(200,*)
	write(200,*) "spin 2 cbm", real(cbm2)
	write(200,*) "spin 2 vbm", real(vbm2)
	write(200,*)
	write(200,*) "spin 2 ponto k cbm:",real(kpt(cbmkpt2,1)),real(kpt(cbmkpt2,2)),real(kpt(cbmkpt2,3))
	write(200,*) "spin 2 ponto k vbm:",real(kpt(vbmkpt2,1)),real(kpt(vbmkpt2,2)),real(kpt(vbmkpt2,3))
	write(200,*)
	write(200,*) "spin 2 ponto k gap direto:",real(kpt(gapdirkpt2,1)),real(kpt(gapdirkpt2,2)),real(kpt(gapdirkpt2,3))
	write(200,*) "spin 2 gp dir cb:", real(energy(gapdirkpt2,nocp2(gapdirkpt2),2))
	write(200,*) "spin 2 gp dir vb:",real(energy(gapdirkpt2,nocp2(gapdirkpt2)-1,2))

	write(200,*) "                  "

	!definindo gap e gap dir completo

	if ( cbm .gt. cbm2) then

		cbmf= cbm2
	else

		cbmf= cbm
	end if

	if ( vbm .gt. vbm2) then

		vbmf= vbm
	else

		vbmf= vbm2
	end if

	write(200,*) "gap indireto", real(cbmf-vbmf)

!	if ( gapdirkpt2 .eq. gapdirkpt) then

!		if (real(energy(gapdirkpt2,nocp2(gapdirkpt2),2)) .gt. real(energy(gapdirkpt,nocp(gapdirkpt),1))) then

!			cbmdir= energy(gapdirkpt,nocp(gapdirkpt),1)
!		else
!			cbmdir= energy(gapdirkpt2,nocp2(gapdirkpt2),2)
!		end if

!		if (real(energy(gapdirkpt2,nocp2(gapdirkpt2)-1,2)) .gt. real(energy(gapdirkpt,nocp(gapdirkpt)-1,1))) then

!			vbmdir= energy(gapdirkpt2,nocp2(gapdirkpt2)-1,2)
!		else
!			vbmdir= energy(gapdirkpt2,nocp2(gapdirkpt2)-1,1)
!		end if

!			write(200,*) "gap direto", real(cbmdir-vbmdir)

!	else

!		if (gapdir .gt. gapdir2) then

!		 	write(200,*) "gap direto s2", real(gapdir2)
!		else
!			write(200,*) "gap direto s1", real(gapdir)
!		end if

!	end if

	gapdirf=500
	gapdirkptf=1

	do i=1,nkpts

	if (energy(i,nocp2(i),2) .gt. energy(i,nocp(i),1)) then

		cbmdir = energy(i,nocp(i),1)

	else 
		cbmdir = energy(i,nocp2(i),2)
	end if

	if (energy(i,nocp2(i)-1,2) .gt. energy(i,nocp(i)-1,1)) then

		vbmdir = energy(i,nocp2(i)-1,2)

	else 
		vbmdir = energy(i,nocp(i)-1,1)
	end if

	gapfaux= cbmdir-vbmdir

	if (gapfaux .lt. gapdirf) then

		gapdirf=gapfaux
		gapdirkptf=i

	else

	end if

	end do

	write(200,*) "gap direto",real(gapdirf)
	write(200,*) "          "
	write(200,*) "ponto k gap direto:",real(kpt(gapdirkptf,1)),real(kpt(gapdirkptf,2)),real(kpt(gapdirkptf,3))

	go to 120

140	write(200,*) "sistema metalico"


120	continue

	deallocate(kpt,nocp,energy,nocp2)
	


	close(100)
	close(200)


end program main
