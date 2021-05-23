!programa para ler a saida das bandas do vasp e gerar o dat para plotar as bandas

!gfortran vasp_bands_kpath_scs.f90 -o vbands_k_scs.x
! vbands_k_scs.x npath nocp scs >bands.dat

program main

	implicit none
	
	integer :: nbands,nkpoints
	double precision :: flag
	double precision,allocatable,dimension(:,:,:) :: energies

	integer :: nocp
	double precision :: scs

	integer :: i,j,k,l,erro

	integer :: nptpath,ntrechos,nptrechos
	double precision,allocatable,dimension(:,:) :: kpt
	double precision,dimension(3,3) :: latvec
	double precision,allocatable,dimension(:) :: kpathp

	!variaveis getarg

  	character*30 arg
  	integer ios

	! INPUT : lendo os parametros do modelo de tight-binding
	OPEN(UNIT=100, FILE= "EIGENVAL",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de entrada 1"
	OPEN(UNIT=101, FILE= "KPOINTS",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de entrada 2"
	OPEN(UNIT=102, FILE= "POSCAR",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de entrada 3"

      CALL getarg(1, arg)
		read(arg,*,iostat=ios) nptpath

      CALL getarg(2, arg)
		read(arg,*,iostat=ios) nocp

      CALL getarg(3, arg)
		read(arg,*,iostat=ios) scs

	ntrechos=nptpath/2

 	!lendo EIGENVAL
	do i=1,5
		read(100,*)
	end do

	read(100,*) flag,nkpoints,nbands

	allocate(energies(nkpoints,nbands,2))

	do i=1,nkpoints

		read(100,*) flag,flag,flag,flag

		do j=1,nbands

			energies(i,j,1) = dble(i)


			read(100,*) flag,energies(i,j,2),flag

			

		end do

		

	end do

	!termino leitura EIGENVAL

	!inicio leitura KPOINTS
	read(101,*)
	read(101,*) nptrechos
	read(101,*)
	read(101,*)

	allocate(kpt(nptpath,3))
	
	do i=1,nptpath

	 read(101,*) kpt(i,1),kpt(i,2),kpt(i,3)

	end do

	!termino leitura KPOINTS

	!inicio leitura POSCAR

	read(102,*)
	read(102,*)
	read(102,*) latvec(1,1),latvec(1,2),latvec(1,3)
	read(102,*) latvec(2,1),latvec(2,2),latvec(2,3)
	read(102,*) latvec(3,1),latvec(3,2),latvec(3,3)


	!termino leitura POSCAR

	allocate(kpathp((nptpath/2)*nptrechos))

	call kpath2(latvec(1,:),latvec(2,:),latvec(3,:),nptpath,kpt,nptrechos,kpathp)


	do j=1,nbands

		do i=1,nkpoints

			if (j .gt. nocp) then

			write(*,*) kpathp(i),energies(i,j,2)+scs

			else

			write(*,*) kpathp(i),energies(i,j,2)

			end if

		end do

			write(*,*)

	end do



	deallocate(energies)
	deallocate(kpt)

	close(100)
	close(101)
	close(102)


end program main

!subrotinas auxiliares

subroutine kpath2(rlat1,rlat2,rlat3,nks,ks,npts,kpt)

	implicit none

	integer :: i,j
	integer :: nks,npts
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(3) :: blat1,blat2,blat3

	double precision,dimension((nks/2)*npts) :: kpt

	double precision,dimension(nks,3) :: ks
	double precision,dimension(nks,3) :: ksaux

	double precision,dimension(nks/2) :: kdis !distancia entre os pontos k do caminho
	double precision :: kdistot !soma do caminho todo
	double precision :: kdisaux


	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	do i=1,nks
	
		ksaux(i,1) = ks(i,1)*blat1(1)+ks(i,2)*blat2(1)+ks(i,3)*blat3(1)
		ksaux(i,2) = ks(i,1)*blat1(2)+ks(i,2)*blat2(2)+ks(i,3)*blat3(2)
		ksaux(i,3) = ks(i,1)*blat1(3)+ks(i,2)*blat2(3)+ks(i,3)*blat3(3)
	
	end do



	kdistot= 0.0
	do i=1,nks/2

		call distvec(ksaux(2*i,:),ksaux(2*i-1,:),kdis(i))

		!kdis(i)=sqrt((ksaux(2*i,1)-ksaux(2*i-1,1))**2+(ksaux(2*i,2)-ksaux(2*i-1,2))**2+(ksaux(2*i,3)-ksaux(2*i-1,3))**2)
		kdistot=kdistot+kdis(i)

	end do

	kdis=kdis/kdistot


	kdisaux=0.0

	 do j=1,nks/2
		
	  do i=1,npts

		!kpt((j-1)*npts+i) = (kdisaux)+(kdis(j))*dble((i)/(npts))
		kpt((j-1)*npts+i) = (kdisaux)+(kdis(j))*(i-1.)/(npts-1.)
		
		
		

	  end do

		kdisaux=kdisaux+kdis(j)
	 end do




end subroutine kpath2

subroutine distvec(vec1,vec2,distvecout) !subroutina para calcular a distância entre dois vetores

	implicit none

	double precision,dimension(3) :: vec1,vec2
	double precision :: distvecout

	distvecout=sqrt((vec1(1)-vec2(1))**2+(vec1(2)-vec2(2))**2+(vec1(3)-vec2(3))**2)


end subroutine distvec

subroutine recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3) !calcula os vetores da rede recíproca, a partir dos vetores da rede real

	implicit none
	
	double precision,parameter :: pi=acos(-1.)
	
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(3) :: blat1,blat2,blat3

	double precision,dimension(3) :: v23,v31,v12
	double precision :: vol

	call prodvec(rlat2,rlat3,v23)
	call prodvec(rlat3,rlat1,v31)
	call prodvec(rlat1,rlat2,v12)

	vol= (rlat1(1)*v23(1))+(rlat1(2)*v23(2))+(rlat1(3)*v23(3))

	blat1 = ((2.0*pi)/vol)*v23
	blat2 = ((2.0*pi)/vol)*v31
	blat3 = ((2.0*pi)/vol)*v12


end subroutine recvec

subroutine prodvec(v1,v2,vx) !calcula o vetor oriundo do produto vetorial entre dois vetores v1 X v2

	implicit none
	
	double precision,dimension(3) :: v1,v2,vx

	vx(1) = (v1(2)*v2(3))-(v1(3)*v2(2))
	vx(2) = (v1(3)*v2(1))-(v1(1)*v2(3))
	vx(3) = (v1(1)*v2(2))-(v1(2)*v2(1))


end subroutine prodvec 
