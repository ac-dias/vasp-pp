!programa para achar o bandgap (semi-metais - materiais com gap pequeno)
!gfortran vasp_gap_vbm-finder-sm.f90 -o gapfsm.x

program main

	implicit none

	integer :: nbands,nkpts
	integer :: erro

	integer :: i,j

	double precision :: flag

	double precision,allocatable,dimension(:,:) :: kpt
	integer,allocatable,dimension(:) :: nocp

	double precision,allocatable,dimension(:,:,:) :: energy

	double precision :: gapdir,cbm,vbm
	integer :: cbmkpt,vbmkpt,gapdirkpt
	double precision :: gaux,vaux,caux

	double precision,parameter :: ocptol=0.05
	double precision,parameter :: ocptol2=0.02

	!variaveis getarg

  	character*30 arg
  	integer ios

	! INPUT : lendo os parametros do modelo de tight-binding
	OPEN(UNIT=100, FILE= "EIGENVAL",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de entrada 1"

	! OUTPUT : lendo os parametros do modelo de tight-binding
	OPEN(UNIT=200, FILE= "gap_vasp.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida 1"

     ! CALL getarg(1, arg)
	!	read(arg,*,iostat=ios) nbands
       !CALL getarg(2, arg)
	!	read(arg,*,iostat=ios) nkpts
 


	!nbands=20
	!nkpts=27



	do i=1,5
		read(100,*)
	end do

	read(100,*) flag,nkpts,nbands

	allocate(kpt(nkpts,3),nocp(nkpts),energy(nkpts,nbands,2))

	!lendo os valores de energia nos pontos k

	do i=1,nkpts

		read(100,*) kpt(i,1),kpt(i,2),kpt(i,3),flag

		do j=1,nbands

			read(100,*) flag,energy(i,j,1),energy(i,j,2)


		end do

		

	end do

	!termino da leitura

	!identificando a primeira banda desocupada de cada ponto k

	do i=1,nkpts


		do j=1,nbands

			if (energy(i,j,2) .lt. 1.0-ocptol2) then

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

	!do i=1,nkpts

	!	if ( (energy(i,nocp(i),2) .gt. 0.0) .and. (energy(i,nocp(i),2) .lt. 1.0-ocptol)) then

	!		go to 70

	!	else 

	!		continue

	!	end if

	!end do
		
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


	write(200,*) "gap direto:",real(gapdir)
	write(200,*) "gap indireto", real(cbm-vbm)
	write(200,*)
	write(200,*) "cbm", real(cbm)
	write(200,*) "vbm", real(vbm)
	write(200,*)
	write(200,*) "ponto k cbm:",real(kpt(cbmkpt,1)),real(kpt(cbmkpt,2)),real(kpt(cbmkpt,3))
	write(200,*) "ponto k vbm:",real(kpt(vbmkpt,1)),real(kpt(vbmkpt,2)),real(kpt(vbmkpt,3))
	write(200,*)
	write(200,*) "ponto k gap direto:",real(kpt(gapdirkpt,1)),real(kpt(gapdirkpt,2)),real(kpt(gapdirkpt,3))

	go to 60

70	write(200,*) "sistema metalico"


60	continue	

	deallocate(kpt,nocp,energy)
	


	close(100)
	close(200)


end program main
