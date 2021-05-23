      PROGRAM VASP_DOS_ANALYSES

      IMPLICIT NONE
      INTEGER, PARAMETER :: q=SELECTED_REAL_KIND(15,307)
      CHARACTER(LEN=80) :: FN1
      CHARACTER(LEN=80) :: FN2
      INTEGER, PARAMETER :: OUTUNIT1=44
      INTEGER, PARAMETER :: OUTUNIT2=45
      INTEGER :: I, J, JJ, K, K0, K1, K2 
      INTEGER :: NATOMS, NEDOS, NPINT, NSPIN, LORBIT, NSTATES
      INTEGER :: NATOMI, NATOMF, NOA0
      INTEGER, ALLOCATABLE, DIMENSION(:) :: AUXK
      REAL(q), PARAMETER :: ZERO = 0.0D0
      REAL(q), PARAMETER :: SIGMAZERO = 0.010D0
      REAL(q) :: ENERGY_MAX, ENERGY_MIN
      REAL(q) :: ENERGY_MAXfermi, ENERGY_MINfermi
      REAL(q) :: ENERGY_MAXint, ENERGY_MINint
      REAL(q) :: EFERMI, SUM1, SUME1, KEY0, KEY1, DELTAE
      REAL(q) :: SIGMAIN, ENERGYJ, STEPENERGY, TEMP, XXX 
      REAL(q), ALLOCATABLE, DIMENSION(:) :: ENERGY, AUX1, AUXE1
      REAL(q), ALLOCATABLE, DIMENSION(:) :: KRONECKER, GAUSSOUT 
      REAL(q), ALLOCATABLE, DIMENSION(:) :: COGDOS, COGLDOSSUM
      REAL(q), ALLOCATABLE, DIMENSION(:,:) :: DOS, ENEDOS, DOSSUM
      REAL(q), ALLOCATABLE, DIMENSION(:,:) :: GAUSSDOS 
      REAL(q), ALLOCATABLE, DIMENSION(:,:) :: LDOSTOT, LDOSINT, LDOSSUM
      REAL(q), ALLOCATABLE, DIMENSION(:,:) :: GAUSSLDOSSUM
      REAL(q), ALLOCATABLE, DIMENSION(:,:) :: COGLDOS, ENELDOSSUM 
      REAL(q), ALLOCATABLE, DIMENSION(:,:,:) :: LDOS, ENELDOS
      REAL(q), ALLOCATABLE, DIMENSION(:,:,:) :: GAUSSLDOS 

      ! Open all files 
      OPEN(10,FILE='DOSCAR',STATUS='OLD')
      OPEN(20,FILE='DOS_0.dat')
      OPEN(21,FILE='DOS_0g.dat')
      OPEN(30,FILE='DOSOUT')

      ! The present version supports only LORBIT = 10 
      ! In principle, it can be extended for different LORBIT very easily 
      WRITE(*,*) 'Provide lorbit - The same as in INCAR'
      READ(*,*) LORBIT
      WRITE(*,*) 
      IF (LORBIT < 10 .OR. LORBIT > 10) THEN
         WRITE(*,*) 'ONLY LORBIT = 10 IS SUPPORTED BY THIS PROGRAM'
         STOP
      ENDIF

      ! Read the NSPIN as provided in INCAR file. 
      WRITE(*,*) 'Provide spin (1 or 2) - The same as in INCAR'
      READ(*,*) NSPIN
      WRITE(*,*)
      IF (NSPIN /= 1 .AND. NSPIN /= 2) THEN
         WRITE(*,*) 'WRONG SPIN'
         STOP
      ENDIF

      ! s, p --> states = 2
      ! s, p, d --> states = 3
      ! s, p, d, f --> states = 4
      ! provide the highest angular momentum as input (s=0, p=1, d=2, f=3)
      ! It will define the size of the array.
      ! calculate the correct number of states
      ! for spin=2
      ! s_up, s_dn, p_up, p_dn, d_up, d_dn, f_up, f_dn  tot_up tot_dn
      ! 1     2     3     4     5     6     7     8     9      10
      ! for spin=1
      ! s, p, d, f  tot
      ! 1  2  3  4  5
      WRITE(*,*) 'Provide the highest l-state written in DOSCAR'
      WRITE(*,*) 'You need to look the number of collumns in DOSCAR'
      WRITE(*,*) 'For example, s=0, p=1, d=2, f=3'
      READ(*,*) NSTATES 
      WRITE(*,*) 
      IF (NSPIN == 2) NSTATES = NSPIN*(NSTATES + 1) + 2
      IF (NSPIN == 1) NSTATES = NSPIN*(NSTATES + 1) + 1

      ! Read the number of atoms in DOSCAR
      ! Read maximum and minimum energies, Efermi, NEDOS
      READ(10,*) NATOMS 
      READ(10,*)
      READ(10,*)
      READ(10,*)
      READ(10,*)
      READ(10,*) ENERGY_MAX, ENERGY_MIN, NEDOS, EFERMI
      ENERGY_MAXfermi = ENERGY_MAX - EFERMI 
      ENERGY_MINfermi =  ENERGY_MIN - EFERMI
      WRITE(*,*) 'Provide minimum energy for integration'
      WRITE(*,*) 'Important for center of gravity calculations'
      WRITE(*,*) 'Suggested value (eV)',ENERGY_MINfermi
      READ(*,*) ENERGY_MINint
      WRITE(*,*)
      IF (ENERGY_MINint <= ENERGY_MINfermi .OR. ENERGY_MINint >= ENERGY_MAXfermi) THEN
         WRITE(*,*) 'WRONG VALUE'
         STOP
      ENDIF  
      WRITE(*,*) 'Provide maximum energy for Integration'
      WRITE(*,*) 'Important for center of gravity calculations'
      WRITE(*,*) 'Suggested value for center of gravity = 0.00 eV'
      READ(*,*) ENERGY_MAXint
      WRITE(*,*)
      IF (ENERGY_MAXint <= ENERGY_MINfermi .OR. ENERGY_MAXint >= ENERGY_MAX) THEN
         WRITE(*,*) 'WRONG VALUE'
         STOP
      ENDIF

      ! Define the smearing parameter. Perform smearing or not.
      WRITE(*,*) 'Provide gaussian smearing (|SIGMA| > 0.010 eV)'
      WRITE(*,*) 'YES (SIGMA > 0), NO (SIGMA < 0)' 
      READ(*,*) SIGMAIN
      IF (ABS(SIGMAIN) < SIGMAZERO) THEN
         WRITE(*,*) 'SIGMA IS VERY SMALL: INCREASE SIGMA'
         STOP
      ENDIF 

      ! allocate all the arrays
      ALLOCATE(ENERGY(1:NEDOS),AUX1(1:NEDOS),AUXE1(1:NEDOS))
      ALLOCATE(GAUSSOUT(1:NEDOS))
      ALLOCATE(KRONECKER(1:NATOMS),COGDOS(1:3*NSPIN))
      ALLOCATE(COGLDOSSUM(1:NSTATES),AUXK(1:NATOMS))
      ALLOCATE(DOS(1:3*NSPIN,1:NEDOS),ENEDOS(1:3*NSPIN,1:NEDOS))
      ALLOCATE(GAUSSDOS(1:3*NSPIN,1:NEDOS))
      ALLOCATE(DOSSUM(1:NSPIN,1:NEDOS),ENELDOSSUM(1:NSTATES,1:NEDOS))
      ALLOCATE(COGLDOS(1:NATOMS,1:NSTATES))
      ALLOCATE(LDOS(1:NATOMS,0:NSTATES,1:NEDOS))
      ALLOCATE(ENELDOS(1:NATOMS,1:NSTATES,1:NEDOS))
      ALLOCATE(GAUSSLDOS(1:NATOMS,0:NSTATES,1:NEDOS))
      ALLOCATE(LDOSSUM(0:NSTATES,1:NEDOS))
      ALLOCATE(GAUSSLDOSSUM(0:NSTATES,1:NEDOS))

      ! Define the list of LDOS to sum up.
      ! Sum of a sequence of sum separated atoms
      WRITE(*,*) 'Do you want to sum few LDOS: YES (N > 0) NO (N < 0)'
      READ(*,*) KEY0
      IF (KEY0 > ZERO) THEN
         ! set up the initial values for delta kronecker (0.0D0 for every atom)
         DO I=1,NATOMS,1
            KRONECKER(I) = 0.0D0
         ENDDO 
         WRITE(*,*) 'Sum a sequence : (N > 0)'
         WRITE(*,*) 'Sum separated atoms: (N < 0)'
         READ(*,*) KEY1
         IF (KEY1 > ZERO) THEN
            WRITE(*,*) 'Provide first and last atom at the same line'
            READ(*,*) NATOMI
            READ(*,*) NATOMF
            !set up delta of kronecker 1.0D0 only for specific atoms 
            DO I=NATOMI,NATOMF,1
               KRONECKER(I) = 1.0D0
            ENDDO
         ELSE
            WRITE(*,*) 'Provide the number of atoms'
            READ(*,*) NOA0
            WRITE(*,*) 'Provide all atoms in the same line'
            ! read all the atoms for the sum dos 
            READ(*,*) (AUXK(I), I=1,NOA0)
            ! set up delta of kronecker 1.0D0 only for specific atoms 
            DO I=1,NOA0,1
               KRONECKER(AUXK(I)) = 1.0D0
            ENDDO
         ENDIF
      ENDIF

!C    Build up the energy grid based on the ENERGY_MIN and ENERGY_MAX
!C    The energy grid is written only with 3 significative numbers. 
!C    This procedure minimize this problem.

      ENERGY = 0.0_q
      DELTAE = (ENERGY_MAX - ENERGY_MIN)/(NEDOS - 1)
      DO J=1,NEDOS
         ENERGY(J) = ENERGY_MIN + (J-1)*DELTAE
      ENDDO

!C    READ DOS total in DOSCAR 

      NPINT = 0
      DO J=1,NEDOS
         READ(10,*) XXX, (DOS(K,J), K=1,NSPIN), & 
         (DOSSUM(K,J), K=1,NSPIN)
         ! SHIFT THE ENERGY SCALE (EFERMI = 0.0 eV)
         ENERGY(J) = ENERGY(J) - EFERMI
         ! BUILD-UP THE NEW FUNCTION
         IF (ENERGY(J) < ZERO) NPINT = NPINT + 1
      ENDDO
      
!C    READ ALL THE LDOS
      
      DO I=1,NATOMS
         READ(10,*)
         DO J=1,NEDOS
            ! it reads the first column, but does not use it.
            ! the first column is the same as energy(j)
            READ(10,*) (LDOS(I,K,J), K=0,(NSTATES - NSPIN))
         ENDDO
      ENDDO

!C    GENERATE TOTAL LDOS FOR EVERY ATOM

      DO I=1,NATOMS
         DO J=1,NEDOS
            IF (NSPIN == 2) THEN
               DO K=1,NSTATES-3,2
                  ! the total ldos are allocated in the highest components in the array
                  ! allocate up component
                  LDOS(I,NSTATES-1,J) = LDOS(I,NSTATES-1,J) + LDOS(I,K,J)
               ENDDO
               DO K=2,NSTATES-2,2
                  ! allocate down component 
                  LDOS(I,NSTATES,J) = LDOS(I,NSTATES,J) + LDOS(I,K,J)
               ENDDO
            ENDIF
            IF (NSPIN == 1) THEN 
               DO K=1,NSTATES-1,1
                  LDOS(I,NSTATES,J) = LDOS(I,NSTATES,J) + LDOS(I,K,J)
               ENDDO
            ENDIF
         ENDDO
      ENDDO            

!C    SUM UP EVERY TOTAL LDOS TO OBTAIN THE TOTAL LDOS

      IF (NSPIN == 2) THEN
         DO J=1,NEDOS,1
            DO I=1,NATOMS,1
               DOS(3,J) = DOS(3,J) + LDOS(I,NSTATES-1,J)
               DOS(4,J) = DOS(4,J) + LDOS(I,NSTATES,J)
            ENDDO
         ENDDO
      ELSE 
         DO J=1,NEDOS,1
            DO I=1,NATOMS,1
               DOS(2,J) = DOS(2,J) + LDOS(I,NSTATES,J)
            ENDDO
         ENDDO
      ENDIF

!C    Calculate the DOSINT and OUTPUT in DOS0.dat
!C    By definition, the sum of the LDOS is smaller than DOS total
!C    Then, by definition DOSINT is always a positive number 

      IF (NSPIN == 2) THEN
         DO J=1,NEDOS
            DOS(5,J) = DOS(1,J) - DOS(3,J)
            DOS(6,J) = DOS(2,J) - DOS(4,J)
         ENDDO
      ELSE
         DO J=1,NEDOS
            DOS(3,J) = DOS(1,J) - DOS(2,J)
         ENDDO                                
      ENDIF

!C    CROSS-CHECK THE DATA based on the condiction that DOSINT > 0 

      K1 = 0
      IF (NSPIN == 2) THEN
         DO J=1,NEDOS
            IF (DOS(5,J) < ZERO) K1 = K1 + 1
            IF (DOS(6,J) < ZERO) K1 = K1 + 2
         ENDDO
         WRITE(*,*) 'ERROR IN DOSCAR. LDOSINT < 0 ',(K1*0.5_q/NEDOS)*100.0_q,'%' 
      ELSE
         DO J=1,NEDOS
            IF (DOS(3,J) < ZERO) K1 = K1 + 1 
         ENDDO
         WRITE(*,*) 'ERROR IN DOSCAR. LDOSINT < 0',(K1/NEDOS)*100.0_q
      ENDIF

!C    BUILD-UP THE NEW FUNCTION  

      DO K=1,3*NSPIN
         DO J=1,NEDOS
            ENEDOS(K,J) = ENERGY(J)*DOS(K,J)
         ENDDO
      ENDDO

!C    CALCULATE THE CENTER OF GRAVITY FOR THE TOTAL DOS

      DO K=1,3*NSPIN
         DO J=1,NEDOS
            IF (ENERGY(J) >= ENERGY_MINint .AND. & 
                ENERGY(J) <= ENERGY_MAXint) THEN
               AUX1(J) = DOS(K,J)
               AUXE1(J) = ENEDOS(K,J)
            ELSE
               AUX1(J) = 0.0D0
               AUXE1(J) = 0.0D0
            ENDIF
         ENDDO
         CALL INTEGRATION_SIMPSON13(1,NEDOS,ENERGY,AUX1,SUM1)
         CALL INTEGRATION_SIMPSON13(1,NEDOS,ENERGY,AUXE1,SUME1)
         IF (NSPIN == 2) THEN
            IF (K == 1) WRITE(*,*) ABS(SUM1 - DOSSUM(1,NPINT))
            IF (K == 2) WRITE(*,*) ABS(SUM1 - DOSSUM(2,NPINT))
         ELSE
            IF (K == 1) WRITE(*,*) ABS(SUM1 - DOSSUM(1,NPINT))
         ENDIF
         ! center of gravity function for total DOS 
         COGDOS(K) = SUME1/SUM1
      ENDDO
      WRITE(30,600) (COGDOS(K), K=1,3*NSPIN)
 600  FORMAT(2X,'0',6F10.4)

!C    PERFORM THE GAUSSIAN SMEARING OF THE DOS FILE

      IF (SIGMAIN > ZERO) THEN  
         ! ENERGY(J) - array with the energy grid
         ! ENERGYJ - particular energy in the grid
         ! SIGMAIN - gaussian smearing parameter (input)
         ! GAUSSOUT - result from the subroutine
         STEPENERGY = ENERGY(2) - ENERGY(1)
         DO J=1,NEDOS
            ENERGYJ = ENERGY(J)
            CALL GAUSS(ENERGY,NEDOS,ENERGYJ,SIGMAIN,GAUSSOUT)
            DO K=1,3*NSPIN
               DO JJ=1,NEDOS
                  TEMP = DOS(K,J)*STEPENERGY*GAUSSOUT(JJ)
                  GAUSSDOS(K,JJ) = GAUSSDOS(K,JJ) + TEMP
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!C    MULTIPLY THE DOS_DOWN BY -1.0 (HELP TO PLOT)

      IF (NSPIN == 2) THEN
         DO J=1,NEDOS
            DO K=1,3*NSPIN
               IF (K==2.OR.K==4.OR.K==6) THEN
                  DOS(K,J) = -1.0D0*DOS(K,J)
                  IF (SIGMAIN > ZERO) GAUSSDOS(K,J)=-1.0D0*GAUSSDOS(K,J)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

!C    OUTPUT TOTAL DOS, TOTAL LDOS, AND INTERSTITAL DOS

      DO J=1,NEDOS
         WRITE(20,300) ENERGY(J), (DOS(K,J), K=1,3*NSPIN)
         IF (SIGMAIN > ZERO) THEN
            WRITE(21,300) ENERGY(J), (GAUSSDOS(K,J), K=1,3*NSPIN)
         ENDIF
 300     FORMAT(F10.3,6E12.4)
      ENDDO

!C    GENERATE DOS SUM FOR A GIVE NUMBER OF ATOMS AND STOP 

      IF (KEY0 > ZERO) THEN
         ! sum-up every LDOS to yield the total LDOS
         DO J=1,NEDOS,1
            DO K=1,(NSTATES - NSPIN),1
               DO I=1,NATOMS,1
                  LDOSSUM(K,J) = LDOSSUM(K,J) + LDOS(I,K,J)*KRONECKER(I)
               ENDDO
            ENDDO
         ENDDO
         IF (NSPIN == 2) THEN
            DO J=1,NEDOS,1
               DO K=1,NSTATES-3,2
                  DO I=1,NATOMS,1 
                     LDOSSUM(NSTATES-1,J) = LDOSSUM(NSTATES-1,J) + & 
                                            LDOS(I,K,J)*KRONECKER(I)
                  ENDDO
               ENDDO
            ENDDO
            DO J=1,NEDOS,1
               DO K=2,NSTATES-2,2
                  DO I=1,NATOMS
                     LDOSSUM(NSTATES,J) = LDOSSUM(NSTATES,J) + & 
                                          LDOS(I,K,J)*KRONECKER(I)
                  ENDDO
               ENDDO
            ENDDO 
         ENDIF
         IF (NSPIN == 1) THEN 
            DO J=1,NEDOS
               DO K=1,NSTATES-1
                  DO I=1,NATOMS 
                     LDOSSUM(NSTATES,J) = LDOSSUM(NSTATES,J) + & 
                                          LDOS(I,K,J)*KRONECKER(I)
                  ENDDO
               ENDDO
            ENDDO 
         ENDIF
         
         ! buld-up the function E*DOS(E)
         DO K=1,NSTATES 
            DO J=1,NEDOS
               ENELDOSSUM(K,J) = ENERGY(J)*LDOSSUM(K,J)
            ENDDO
         ENDDO
         ! calculate the center of gravity for every state
         DO K=1,NSTATES
            DO J=1,NEDOS
               ! integration is performed only for a defined range.
               IF (ENERGY(J) >= ENERGY_MINint .AND. & 
                   ENERGY(J) <= ENERGY_MAXint) THEN
                  AUX1(J) = LDOSSUM(K,J)
                  AUXE1(J) = ENELDOSSUM(K,J)
               ELSE
                  AUX1(J) = 0.0D0
                  AUXE1(J) = 0.0D0
               ENDIF
            ENDDO
            ! Integration of the two functions.  
            CALL INTEGRATION_SIMPSON13(1,NEDOS,ENERGY,AUX1,SUM1)
            CALL INTEGRATION_SIMPSON13(1,NEDOS,ENERGY,AUXE1,SUME1)
            IF (SUM1 > ZERO) THEN
               COGLDOSSUM(K) = SUME1/SUM1
            ELSE
               COGLDOSSUM(K) = 99
            ENDIF
         ENDDO
         ! write the center of gravity for every state
         K0 = NSTATES
         K1 = NSTATES - NSPIN + 1
         K2 = NSTATES - NSPIN 
         WRITE(30,700) (COGLDOSSUM(K), K=K1,K0), (COGLDOSSUM(K), K=1,K2)
 700     FORMAT('SUM',10F10.4)
         ! Gaussian smearing
         IF (SIGMAIN > ZERO) THEN
            STEPENERGY = ENERGY(2) - ENERGY(1)
            DO J=1,NEDOS
               ENERGYJ = ENERGY(J)
               CALL GAUSS(ENERGY,NEDOS,ENERGYJ,SIGMAIN,GAUSSOUT)
               DO K=1,NSTATES
                  DO JJ=1,NEDOS
                     TEMP = LDOSSUM(K,J)*STEPENERGY*GAUSSOUT(JJ)
                     GAUSSLDOSSUM(K,JJ) = GAUSSLDOSSUM(K,JJ) + TEMP
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         ! add the minus sign for the down components
         IF (NSPIN == 2) THEN
            DO K=1,NSTATES
               DO J=1,NEDOS
                  IF (K==2 .OR. K==4 .OR. K==6 .OR. K==8 .OR. K==10) THEN 
                     LDOSSUM(K,J) = -1.0D0*LDOSSUM(K,J)
                     IF (SIGMAIN > ZERO) THEN
                        GAUSSLDOSSUM(K,J) = -1.0D0*GAUSSLDOSSUM(K,J)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         ! output the total ldos sum.
         OPEN(40,FILE='DOS_sum.dat')
         OPEN(41,FILE='DOS_sumg.dat')
         ! tot_up tot_dn s_up s_dn p_up p_dn d_up d_dn f_up f_dn
         DO J=1,NEDOS
            WRITE(40,100) ENERGY(J), (LDOSSUM(K,J), K=K1,K0), & 
                 (LDOSSUM(K,J), K=1,K2)
            IF (SIGMAIN > ZERO) THEN
               WRITE(41,100) ENERGY(J), & 
               (GAUSSLDOSSUM(K,J), K=K1,K0), (GAUSSLDOSSUM(K,J), K=1,K2)
            ENDIF
         ENDDO 
         STOP 
      ENDIF
  
!C    BUILD-UP NEW FUNCTION FOR CENTER OF GRAVITY LDOS CALCULATION 

      DO I=1,NATOMS
         DO K=1,NSTATES  
            DO J=1,NEDOS 
               ENELDOS(I,K,J) = ENERGY(J)*LDOS(I,K,J)
            ENDDO
         ENDDO
      ENDDO
      
!C    NUMERICAL INTEGRATION TO OBTAIN THE CENTER OF GRAVITY FOR EVERY FUNCTION 

      DO I=1,NATOMS
         DO K=1,NSTATES
            DO J=1,NEDOS
               ! integration is performed only for a defined range.
               IF (ENERGY(J) >= ENERGY_MINint .AND. & 
                   ENERGY(J) <= ENERGY_MAXint) THEN
                  AUX1(J) = LDOS(I,K,J)
                  AUXE1(J) = ENELDOS(I,K,J)
               ELSE
                  AUX1(J) = 0.0D0
                  AUXE1(J) = 0.0D0
               ENDIF
            ENDDO
            ! Integration of the two functions.  
            CALL INTEGRATION_SIMPSON13(1,NEDOS,ENERGY,AUX1,SUM1)
            CALL INTEGRATION_SIMPSON13(1,NEDOS,ENERGY,AUXE1,SUME1)
            IF (SUM1 > ZERO) THEN 
               COGLDOS(I,K) = SUME1/SUM1
            ELSE
               COGLDOS(I,K) = 99 
            ENDIF
         ENDDO
      ENDDO

!C    OUTPUT ALL THE CENTER OF GRAVITY FOR EVERY ATOM AND STATE
          
      !ENE  tot_up tot_dn s_up s_dn p_up p_dn d_up d_dn f_up f_dn
      DO I=1,NATOMS
         K0 = NSTATES
         K1 = NSTATES - NSPIN + 1
         K2 = NSTATES - NSPIN
         WRITE(30,500) I, (COGLDOS(I,K), K=K1,K0), &
                          (COGLDOS(I,K), K=1,K2)
 500     FORMAT(I3,10F10.4)
      ENDDO

!C    PERFORM THE GAUSSIAN SMEARING OF THE LDOS FILE

      IF (SIGMAIN > ZERO) THEN
         ! ENERGY(J) - array with the energy grid
         ! ENERGYJ - particular energy in the grid
         ! SIGMAIN - gaussian smearing parameter (input)
         ! GAUSSOUT - result from the subroutine
         STEPENERGY = ENERGY(2) - ENERGY(1)
         DO I=1,NATOMS
            DO J=1,NEDOS
               ENERGYJ = ENERGY(J)
               CALL GAUSS(ENERGY,NEDOS,ENERGYJ,SIGMAIN,GAUSSOUT)
               DO K=1,NSTATES
                  DO JJ=1,NEDOS
                     TEMP = LDOS(I,K,J)*STEPENERGY*GAUSSOUT(JJ)
                     GAUSSLDOS(I,K,JJ) = GAUSSLDOS(I,K,JJ) + TEMP
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!C    MULTIPLY THE LDOS_DOWN BY -1.0 (HELP TO PLOT)

      IF (NSPIN == 2) THEN
         DO I=1,NATOMS
            DO K=1,NSTATES
               DO J=1,NEDOS
                  IF (K==2 .OR. K==4 .OR. K==6 .OR. K==8 .OR. K==10) THEN 
                     LDOS(I,K,J) = -1.0D0*LDOS(I,K,J)
                     IF (SIGMAIN > ZERO) THEN 
                        GAUSSLDOS(I,K,J) = -1.0D0*GAUSSLDOS(I,K,J)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF 

!C    SEPARATE ALL THE LDOS IN DIFFERENT FILES NUMBERED BY THE ATOMS NUMBERS

      DO I=1,NATOMS
         ! BUILD FILNEMA -- DOS_xxx.DAT
         WRITE(FN1,FMT='(A,I0,A)') 'DOS_',I,'.dat'
         ! OPEN IT WITH A FIXED UNIT NUMBER
         OPEN(UNIT=OUTUNIT1,FILE=FN1,FORM='FORMATTED')
         ! general output for spin=2 and spin=1
         ! energy, tot_up, top_dn, s_up, s_dn, p_up, p_dn, d_up, d_dn, f_up, f_dn 
         DO J=1,NEDOS
            K0 = NSTATES
            K1 = NSTATES - NSPIN + 1
            K2 = NSTATES - NSPIN 
            WRITE(OUTUNIT1,100) ENERGY(J), (LDOS(I,K,J), K=K1,K0), & 
                    (LDOS(I,K,J), K=1,K2)
 100        FORMAT(F10.3,10E12.4)
         ENDDO
         CLOSE(OUTUNIT1)
      ENDDO

      IF (SIGMAIN > ZERO) THEN
         DO I=1,NATOMS
            ! BUILD FILNEMA -- DOS_xxx.DAT
            WRITE(FN2,FMT='(A,I0,A)') 'DOS_',I,'g.dat'
            ! OPEN IT WITH A FIXED UNIT NUMBER
            OPEN(UNIT=OUTUNIT2,FILE=FN2,FORM='FORMATTED')
            ! general output for spin=2 and spin=1
            ! energy, tot_up, top_dn, s_up, s_dn, p_up, p_dn, d_up, d_dn, f_up, f_dn 
            DO J=1,NEDOS
               K0 = NSTATES
               K1 = NSTATES - NSPIN + 1
               K2 = NSTATES - NSPIN
               WRITE(OUTUNIT2,100) ENERGY(J), & 
                (GAUSSLDOS(I,K,J), K=K1,K0), (GAUSSLDOS(I,K,J), K=1,K2)
            ENDDO
            CLOSE(OUTUNIT2)
         ENDDO
      ENDIF

      DEALLOCATE(ENERGY,AUX1,AUXE1)
      DEALLOCATE(GAUSSOUT)
      DEALLOCATE(KRONECKER,COGDOS)
      DEALLOCATE(COGLDOSSUM,AUXK)
      DEALLOCATE(DOS,ENEDOS)
      DEALLOCATE(GAUSSDOS)
      DEALLOCATE(DOSSUM,ENELDOSSUM)
      DEALLOCATE(COGLDOS)
      DEALLOCATE(LDOS)
      DEALLOCATE(ENELDOS)
      DEALLOCATE(GAUSSLDOS)
      DEALLOCATE(LDOSSUM)
      DEALLOCATE(GAUSSLDOSSUM)

      ! stop the main program.
      STOP

      ! This program include the following subroutines below
      CONTAINS

      SUBROUTINE INTEGRATION_SIMPSON13(A,B,X,F,S)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: A, B
      REAL(q), DIMENSION(:), INTENT(IN) :: X, F
      REAL(q), INTENT(OUT) :: S
      INTEGER :: I, N
      REAL(q) :: S0, S1, S2, H
      N = B - A + 1
      H = (X(B) - X(A))/DBLE(N)
      S0 = 0.0D0
      S1 = 0.0D1
      S2 = 0.0D2
      DO I=A+1,B-1,2
         S1 = S1 + F(I-1)
         S0 = S0 + F(I)
         S2 = S2 + F(I+1) 
      ENDDO
         S = H*(S1 + 4.0D0*S0 + S2)/3.0D0
         ! if n is even, add last slice separately
         IF (MOD(N,2) .EQ. 0) THEN
            S = S + H*(F(N-2) + 4.0D0*F(N-1) + F(N))/3.0D0
         ENDIF
      END SUBROUTINE INTEGRATION_SIMPSON13

      SUBROUTINE GAUSS(Y,N,Y0,SIG,GAUSSOUT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL(q), INTENT(IN) :: Y0, SIG
      REAL(q), DIMENSION(:), INTENT(IN) :: Y
      REAL(q), DIMENSION(:), INTENT(OUT) :: GAUSSOUT
      INTEGER :: J
      REAL(q) :: PI, SIGMA, RNORM, DENO 
      REAL :: TEMP1
     
      PI = 4.0D0*DATAN(1.0D0)
      SIGMA = SIG
      RNORM = 1.0D0/(SIGMA*DSQRT(2.0D0*PI))
      DENO = 2.0D0*SIGMA*SIGMA
      DO J=1,N
         TEMP1 = -1.0D0*(Y(J) - Y0)*(Y(J) - Y0)/DENO
         GAUSSOUT(J) = RNORM*EXP(TEMP1)
      ENDDO
      END SUBROUTINE GAUSS 

      END PROGRAM VASP_DOS_ANALYSES
