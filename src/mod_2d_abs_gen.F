!------------------------------------------------------------------------------------
!
!      FILE mod_2d_abs_gen.F
!
!      This file is part of the FUNWAVE-TVD program under the Simplified BSD license
!
!-------------------------------------------------------------------------------------
! 
!    Copyright (c) 2016, FUNWAVE Development Team
!
!    (See http://www.udel.edu/kirby/programs/funwave/funwave.html
!     for Development Team membership)
!
!    All rights reserved.
!
!    FUNWAVE_TVD is free software: you can redistribute it and/or modify
!    it under the terms of the Simplified BSD License as released by
!    the Berkeley Software Distribution (BSD).
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions are met:
!
!    1. Redistributions of source code must retain the above copyright notice, this
!       list of conditions and the following disclaimer.
!    2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!  
!    The views and conclusions contained in the software and documentation are those
!    of the authors and should not be interpreted as representing official policies,
!    either expressed or implied, of the FreeBSD Project.
!  
!-------------------------------------------------------------------------------------
!
!  TIDE_MODULE is a module to add tide/surge boundary conditions into wave simulation    
!
!  HISTORY :
!    03/25/2025 Fengyan Shi
!-------------------------------------------------------------------------------------

# if defined (MAP2D_ABS_GEN)

MODULE ABS_GEN_2D_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Mloc1,Nloc1,Nghost,Ibeg,Iend,Jbeg,Jend,DX,DY, &
                    H,ETA,U,V,HU,HV,MinDepth,MASK,DT,Depth,tmp4preview, &
                    MASK9,PERIODIC, &
        ! to remove the following sponges
                   Sponge_west_width,&
                   Sponge_east_width,&
                   Sponge_south_width,&
                   Sponge_North_width,&
                   Mglob,Nglob, &
                   TIME,OUT_MAP2D_ETA,OUT_MAP2D_U,OUT_MAP2D_V,RESULT_FOLDER, &
                   WAVEMAKER
  USE INPUT_READ
  USE GLOBAL,ONLY : Cm_eta,Sm_eta,Cm_u,Sm_u,Cm_v,Sm_v,&
                    NumFreq,NumDir,&
                    Per_Ser,Theta_Ser,Phase_LEFT, Amp_Ser, &
                    Phase_Ser,Theta_2D, &
                    Dep_Ser,Segma_Ser,Wave_Number_Ser,Stokes_Drift_Ser

#if defined (PARALLEL)
  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY
  USE MPI
# endif
  IMPLICIT NONE
  SAVE

  CHARACTER(LEN=80) MappingDataFileName
  CHARACTER(LEN=80) TMP_NAME1
  REAL(SP),DIMENSION(:,:),ALLOCATABLE :: SPONGE_2D
  INTEGER :: Ifile,Iwidth_2D_Mapping
  REAL(SP) :: Time_mapping_1,Time_mapping_2
  REAL(SP),DIMENSION(:,:),ALLOCATABLE :: U_2D_INPUT_1,V_2D_INPUT_1,ETA_2D_INPUT_1, &
              U_2D_INPUT_2,V_2D_INPUT_2,ETA_2D_INPUT_2, &
              U_2D_INPUT,V_2D_INPUT,ETA_2D_INPUT
  REAL(SP),DIMENSION(:,:),ALLOCATABLE :: U_2D_MAP,V_2D_MAP,ETA_2D_MAP
                                  
  LOGICAL :: MAP2D_GEN_ABS
  LOGICAL :: OUT_SPONGE_MAP
  LOGICAL :: BC_WEST_NEST,BC_EAST_NEST,BC_SOUTH_NEST,BC_NORTH_NEST

! mapping variables
  INTEGER :: M_DATA,N_DATA
  REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Sc,S1,S2,S3
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: nx1,ny1,nx2,ny2,nx3,ny3
  REAL(SP),DIMENSION(:,:),ALLOCATABLE :: x_grid1,y_grid1,x_grid2,y_grid2

CONTAINS
  
SUBROUTINE MAPPING_INITIAL
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL,INPUT_FILE_NAME
# if defined (PARALLEL)
  USE GLOBAL,ONLY : iista,jjsta   
# endif
                    
  USE INPUT_READ
  IMPLICIT NONE

  CHARACTER(LEN=80)::FILE_NAME=' '
  INTEGER :: Ifile,ierr


! read  from input.txt
      FILE_NAME=INPUT_FILE_NAME

 ! generating and absorbing

      CALL READ_LOGICAL(MAP2D_GEN_ABS,FILE_NAME,'MAP2D_GEN_ABS',ierr)
      IF(ierr == 1)THEN
       MAP2D_GEN_ABS = .FALSE. 
      ELSE
        IF(MAP2D_GEN_ABS)THEN

# if defined (PARALLEL) 
        if (myid.eq.0)then
         WRITE(3,'(A80)')'MAP2D_GEN_ABS is TRUE'
         WRITE(*,'(A80)')'MAP2D_GEN_ABS is TRUE'
        endif
# else
        WRITE(3,'(A80)')'MAP2D_GEN_ABS is TRUE'
        WRITE(*,'(A80)')'MAP2D_GEN_ABS is TRUE'
# endif 
       ENDIF ! true
      ENDIF

! --- nest boundary true or false

      CALL READ_LOGICAL(BC_WEST_NEST,FILE_NAME,'BC_WEST_NEST',ierr)
      IF(ierr == 1)THEN
       BC_WEST_NEST = .FALSE. 
      ELSE
        IF(BC_WEST_NEST)THEN

# if defined (PARALLEL) 
        if (myid.eq.0)then
         WRITE(3,'(A80)')'BC_WEST_NEST is TRUE'
         WRITE(*,'(A80)')'BC_WEST_NEST is TRUE'
        endif
# else
        WRITE(3,'(A80)')'BC_WEST_NEST is TRUE'
        WRITE(*,'(A80)')'BC_WEST_NEST is TRUE'
# endif 
       ENDIF ! true
      ENDIF

      CALL READ_LOGICAL(BC_EAST_NEST,FILE_NAME,'BC_EAST_NEST',ierr)
      IF(ierr == 1)THEN
       BC_EAST_NEST = .FALSE. 
      ELSE
        IF(BC_EAST_NEST)THEN

# if defined (PARALLEL) 
        if (myid.eq.0)then
         WRITE(3,'(A80)')'BC_EAST_NEST is TRUE'
         WRITE(*,'(A80)')'BC_EAST_NEST is TRUE'
        endif
# else
        WRITE(3,'(A80)')'BC_EAST_NEST is TRUE'
        WRITE(*,'(A80)')'BC_EAST_NEST is TRUE'
# endif 
       ENDIF ! true
      ENDIF

      CALL READ_LOGICAL(BC_SOUTH_NEST,FILE_NAME,'BC_SOUTH_NEST',ierr)
      IF(ierr == 1)THEN
       BC_SOUTH_NEST = .FALSE. 
      ELSE
        IF(BC_SOUTH_NEST)THEN

# if defined (PARALLEL) 
        if (myid.eq.0)then
         WRITE(3,'(A80)')'BC_SOUTH_NEST is TRUE'
         WRITE(*,'(A80)')'BC_SOUTH_NEST is TRUE'
        endif
# else
        WRITE(3,'(A80)')'BC_SOUTH_NEST is TRUE'
        WRITE(*,'(A80)')'BC_SOUTH_NEST is TRUE'
# endif 
       ENDIF ! true
      ENDIF

      CALL READ_LOGICAL(BC_NORTH_NEST,FILE_NAME,'BC_NORTH_NEST',ierr)
      IF(ierr == 1)THEN
       BC_NORTH_NEST = .FALSE. 
      ELSE
        IF(BC_NORTH_NEST)THEN

# if defined (PARALLEL) 
        if (myid.eq.0)then
         WRITE(3,'(A80)')'BC_NORTH_NEST is TRUE'
         WRITE(*,'(A80)')'BC_NORTH_NEST is TRUE'
        endif
# else
        WRITE(3,'(A80)')'BC_NORTH_NEST is TRUE'
        WRITE(*,'(A80)')'BC_NORTH_NEST is TRUE'
# endif 
       ENDIF ! true
      ENDIF

! ---

   IF(MAP2D_GEN_ABS)THEN

! width of sponge

      CALL READ_INTEGER(Iwidth_2D_Mapping,FILE_NAME,'Iwidth_2D_Mapping',ierr)

      IF(ierr==1)THEN
        Iwidth_2D_Mapping = 20
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Iwidth_2D_Mapping Default:  20 points'
         WRITE(3,'(A40)')'Iwidth_2D_Mapping Default:  20 points'
      endif
# else
         WRITE(*,'(A40)')'Iwidth_2D_Mapping Default:  20 points'
         WRITE(3,'(A40)')'Iwidth_2D_Mapping Default:  20 points'
# endif
       ENDIF     

      CALL READ_LOGICAL(OUT_MAP2D_ETA,FILE_NAME,'MAP2D_ETA',ierr) 
       IF(OUT_MAP2D_ETA)THEN
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'OUTPUT MAP2D_ETA'
         WRITE(3,'(A40)')'OUTPUT MAP2D_ETA'
      endif
# else
         WRITE(*,'(A40)')'OUTPUT MAP2D_ETA'
         WRITE(3,'(A40)')'OUTPUT MAP2D_ETA'
# endif
       ENDIF
 
      CALL READ_LOGICAL(OUT_MAP2D_U,FILE_NAME,'MAP2D_U',ierr) 
       IF(OUT_MAP2D_U)THEN
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'OUTPUT MAP2D_U'
         WRITE(3,'(A40)')'OUTPUT MAP2D_U'
      endif
# else
         WRITE(*,'(A40)')'OUTPUT MAP2D_U'
         WRITE(3,'(A40)')'OUTPUT MAP2D_U'
# endif
       ENDIF


      CALL READ_LOGICAL(OUT_MAP2D_V,FILE_NAME,'MAP2D_V',ierr) 
       IF(OUT_MAP2D_V)THEN
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'OUTPUT MAP2D_V'
         WRITE(3,'(A40)')'OUTPUT MAP2D_V'
      endif
# else
         WRITE(*,'(A40)')'OUTPUT MAP2D_V'
         WRITE(3,'(A40)')'OUTPUT MAP2D_V'
# endif
       ENDIF

      CALL READ_LOGICAL(OUT_SPONGE_MAP,FILE_NAME,'SPONGE_MAP',ierr) 
       IF(OUT_SPONGE_MAP)THEN
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'OUTPUT SPONGE_MAP'
         WRITE(3,'(A40)')'OUTPUT SPONGE_MAP'
      endif
# else
         WRITE(*,'(A40)')'OUTPUT SPONGE_MAP'
         WRITE(3,'(A40)')'OUTPUT SPONGE_MAP'
# endif
       ENDIF

      ALLOCATE(SPONGE_2D(Mloc,Nloc),U_2D_MAP(Mloc,Nloc), &
               V_2D_MAP(Mloc,Nloc),ETA_2D_MAP(Mloc,Nloc))


! tmp!!!


! 1) ------
      CALL MAP2D_READ_DATA_INIT

! 2) ------
      CALL MAP2D_SPONGE

      CALL REMOVE_ALL_SPONGE

   ENDIF ! end if MAP2D_GEN_ABS

END SUBROUTINE MAPPING_INITIAL


SUBROUTINE MAP2D_SPONGE
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL
# if defined (PARALLEL)
  USE GLOBAL,ONLY : iista,jjsta   
# endif
                    
  IMPLICIT NONE
  REAL(SP) :: ri,lim,Lstart,Lend,R_sp,A_sp
  REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Sponge_w,Sponge_e,Sponge_s,Sponge_n

  ALLOCATE(Sponge_w(Mloc,Nloc),Sponge_e(Mloc,Nloc), &
           Sponge_s(Mloc,Nloc),Sponge_n(Mloc,Nloc))
  
  lim = 1.0_SP
  R_sp = 0.85_SP
  R_sp = 0.95_SP
  A_sp = 10.0_SP

! west

       do j = 1,Nloc
       do i = 1,Mloc
# if defined(PARALLEL)
         ri = R_sp**(50*(i+npx*Mglob/px-1)/(Iwidth_2D_Mapping-1))
# else
         ri = R_sp**(50*(i-1)/(Iwidth_2D_Mapping-1))
# endif
         SPONGE_w(i,j) = max(A_sp**ri,lim)
       enddo
       enddo

! east
       do j = 1,Nloc
       do i = 1,Mloc
# if defined(PARALLEL)
         ri = R_sp**(50*(Mloc-i+(px-npx-1)*Mglob/px)/(Iwidth_2D_Mapping-1))
# else
         ri = R_sp**(50*(Mloc-i)/(Iwidth_2D_Mapping-1))
# endif
         SPONGE_e(i,j) = max(A_sp**ri,lim)
       enddo
       enddo

! south

       do j = 1,Nloc
       do i = 1,Mloc
# if defined(PARALLEL)
         ri = R_sp**(50*(j+npy*Nglob/py-1)/(Iwidth_2D_Mapping-1))
# else
         ri = R_sp**(50*(j-1)/(Iwidth_2D_Mapping-1))
# endif
         SPONGE_s(i,j) = max(A_sp**ri,lim)
       enddo
       enddo

! north
       do j = 1,Nloc
       do i = 1,Mloc
# if defined(PARALLEL)
         ri = R_sp**(50*(Nloc-j+(py-npy-1)*Nglob/py)/(Iwidth_2D_Mapping-1))
# else
         ri = R_sp**(50*(Nloc-j)/(Iwidth_2D_Mapping-1))
# endif
         SPONGE_n(i,j) = max(A_sp**ri,lim)
       enddo
       enddo

! the make values from 0 to 1
      IF(BC_WEST_NEST)THEN
       SPONGE_w = 1.0_SP/SPONGE_w
      ELSE
       SPONGE_w = 1.0_SP
      ENDIF

      IF(BC_EAST_NEST)THEN
       SPONGE_e = 1.0_SP/SPONGE_e
      ELSE
       SPONGE_e = 1.0_SP
      ENDIF

      IF(BC_SOUTH_NEST)THEN
       SPONGE_s = 1.0_SP/SPONGE_s
      ELSE
       SPONGE_s = 1.0_SP
      ENDIF

      IF(BC_NORTH_NEST)THEN
       SPONGE_n = 1.0_SP/SPONGE_n
      ELSE
       SPONGE_n = 1.0_SP
      ENDIF
 
       DO j = 1,Nloc
       DO i = 1,Mloc
         SPONGE_2D(I,J) = MIN(SPONGE_w(I,J),SPONGE_e(I,J),SPONGE_s(I,J),SPONGE_n(I,J))
       ENDDO
       ENDDO

       IF(OUT_SPONGE_MAP)THEN
        TMP_NAME1 = TRIM(RESULT_FOLDER)//'SpongeMap.txt'
        call PutFile(TMP_NAME1,SPONGE_2D)
       ENDIF

       DEALLOCATE(Sponge_w,Sponge_e,Sponge_s,Sponge_n)


END SUBROUTINE  MAP2D_SPONGE

SUBROUTINE REMOVE_ALL_SPONGE
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL,INPUT_FILE_NAME,WaveMaker
# if defined (PARALLEL)
  USE GLOBAL,ONLY : iista,jjsta   
# endif
 
  IMPLICIT NONE

        Sponge_west_width = ZERO
        Sponge_east_width = ZERO
        Sponge_south_width = ZERO
        Sponge_North_width = ZERO

END SUBROUTINE REMOVE_ALL_SPONGE

SUBROUTINE MAP2D_READ_DATA_INIT
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL,INPUT_FILE_NAME,WaveMaker
# if defined (PARALLEL)
  USE GLOBAL,ONLY : iista,jjsta   
# endif
                    
  USE INPUT_READ
  IMPLICIT NONE

  CHARACTER(LEN=80)::FILE_NAME=' '
  INTEGER :: Ifile,ierr
! read  from input.txt
      FILE_NAME=INPUT_FILE_NAME


! read data file name -----------
       CALL READ_STRING(MappingDataFileName,FILE_NAME,'MappingDataFileName',ierr)

        IF(ierr==1)THEN
# if defined (PARALLEL)
        IF(MYID==0)  &
         WRITE(*,*) TRIM(TMP_NAME1), ' No MappingDataFileName, STOP!'
        call MPI_FINALIZE ( ier )
# else
         WRITE(*,*) TRIM(TMP_NAME1), ' No MappingDataFileName, STOP!'        
# endif
         STOP
        ELSE
          TMP_NAME1 = TRIM(MappingDataFileName)
        ENDIF
      
        Ifile = 301
        OPEN(Ifile,FILE=TRIM(TMP_NAME1))   !!! here

      READ(Ifile,'(A80)')  TMP_NAME1 
      READ(Ifile,*)M_DATA,N_DATA

# if defined (PARALLEL)
   IF(MYID==0)THEN
     WRITE(*,*)'DATA dimensions: ', M_DATA,N_DATA
     WRITE(3,*)'DATA dimensions: ', M_DATA,N_DATA
   ENDIF
# else
     WRITE(*,*)'DATA dimensions: ', M_DATA,N_DATA
     WRITE(3,*)'DATA dimensions: ', M_DATA,N_DATA
# endif
      


      ALLOCATE(U_2D_INPUT_1(M_DATA,N_DATA),V_2D_INPUT_1(M_DATA,N_DATA), &
               ETA_2D_INPUT_1(M_DATA,N_DATA), &
               U_2D_INPUT_2(M_DATA,N_DATA),V_2D_INPUT_2(M_DATA,N_DATA), &
               ETA_2D_INPUT_2(M_DATA,N_DATA), &
               U_2D_INPUT(M_DATA,N_DATA),V_2D_INPUT(M_DATA,N_DATA), &
               ETA_2D_INPUT(M_DATA,N_DATA))


      ALLOCATE(Sc(Mloc,Nloc),S1(Mloc,Nloc),S2(Mloc,Nloc),S3(Mloc,Nloc), &
             nx1(Mloc,Nloc),ny1(Mloc,Nloc),nx2(Mloc,Nloc), &
             ny2(Mloc,Nloc),nx3(Mloc,Nloc),ny3(Mloc,Nloc))
      ALLOCATE(x_grid1(M_DATA,N_DATA),y_grid1(M_DATA,N_DATA))
               
      ALLOCATE(x_grid2(Mloc,Nloc),y_grid2(Mloc,Nloc))

      READ(Ifile,'(A80)')  TMP_NAME1  ! title x coordinate

        DO J=1,N_DATA
          READ(Ifile,*) (x_grid1(I,J),I=1,M_DATA)
        ENDDO

# if defined (PARALLEL)
   IF(MYID==0)THEN
     WRITE(*,*)'data X coordinate read over'
     WRITE(3,*)'data X coordinate read over'
   ENDIF
# else
     WRITE(*,*)'data X coordinate read over'
     WRITE(3,*)'data X coordinate read over'
# endif

      READ(Ifile,'(A80)')  TMP_NAME1  ! title y coordinate

        DO J=1,N_DATA
          READ(Ifile,*) (y_grid1(I,J),I=1,M_DATA)
        ENDDO

# if defined (PARALLEL)
   IF(MYID==0)THEN
     WRITE(*,*)'data Y coordinate read over'
     WRITE(3,*)'data Y coordinate read over'
   ENDIF
# else
     WRITE(*,*)'data Y coordinate read over'
     WRITE(3,*)'data Y coordinate read over'
# endif
        
        CALL MAPPING_COEFF     

      DEALLOCATE(x_grid1,y_grid1)
      DEALLOCATE(x_grid2,y_grid2)       

        READ(Ifile,*) Time_mapping_2
        DO J=1,N_DATA
          READ(Ifile,*) (ETA_2D_INPUT_2(I,J),I=1,M_DATA)
        ENDDO
        DO J=1,N_DATA
          READ(Ifile,*) (U_2D_INPUT_2(I,J),I=1,M_DATA)
        ENDDO
        DO J=1,N_DATA
          READ(Ifile,*) (V_2D_INPUT_2(I,J),I=1,M_DATA)
        ENDDO

        Time_mapping_1 = Time_mapping_2
        ETA_2D_INPUT_1 = ETA_2D_INPUT_2
        U_2D_INPUT_1 = U_2D_INPUT_2
        V_2D_INPUT_1 = V_2D_INPUT_2


# if defined (PARALLEL)
   IF(MYID==0)THEN
     WRITE(*,*)'Read mapping data ------'
     WRITE(*,*)'Time = ', Time_mapping_2
     WRITE(3,*)'Read mapping data ------'
     WRITE(3,*)'Time = ', Time_mapping_2
   ENDIF
# else
     WRITE(*,*)'Read mapping data ------'
     WRITE(*,*)'Time = ', Time_mapping_2
     WRITE(3,*)'Read mapping data ------'
     WRITE(3,*)'Time = ', Time_mapping_2
# endif


END SUBROUTINE MAP2D_READ_DATA_INIT

SUBROUTINE MAP2D_READ_DATA
    USE GLOBAL,ONLY : tmp1,tmp2,SMALL,TIME,ZERO,DT
        INTEGER :: Ifile


! west ------------------
        Ifile = 301  

        IF(TIME>Time_mapping_1.AND.TIME>Time_mapping_2) THEN

        Time_mapping_1 = Time_mapping_2
        ETA_2D_INPUT_1 = ETA_2D_INPUT_2
        U_2D_INPUT_1 = U_2D_INPUT_2
        V_2D_INPUT_1 = V_2D_INPUT_2

        DO WHILE (Time_mapping_2.LT.TIME+DT)

        READ(Ifile,*) Time_mapping_2
        DO J=1,N_DATA
          READ(Ifile,*) (ETA_2D_INPUT_2(I,J),I=1,M_DATA)
        ENDDO
        DO J=1,N_DATA
          READ(Ifile,*) (U_2D_INPUT_2(I,J),I=1,M_DATA)
        ENDDO
        DO J=1,N_DATA
          READ(Ifile,*) (V_2D_INPUT_2(I,J),I=1,M_DATA)
        ENDDO

        END DO

# if defined (PARALLEL)
   IF(MYID==0)THEN
     WRITE(*,*)'Read mapping data ------'
     WRITE(*,*)'Time = ', Time_mapping_2
     WRITE(3,*)'Read mapping data ------'
     WRITE(3,*)'Time = ', Time_mapping_2
   ENDIF
# else
     WRITE(*,*)'Read mapping data ------'
     WRITE(*,*)'Time = ', Time_mapping_2
     WRITE(3,*)'Read mapping data ------'
     WRITE(3,*)'Time = ', Time_mapping_2
# endif

     ENDIF ! end read time

! interpolate
    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>Time_mapping_1)THEN
      IF(Time_mapping_1.EQ.Time_mapping_2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(Time_mapping_2-TIME) &
            /MAX(SMALL, ABS(Time_mapping_2-Time_mapping_1))
      tmp1=1.0_SP - tmp2
      ENDIF  ! no more data?
    ENDIF ! time>time_1

! interpolation in time
  
   DO J=1,N_DATA
   DO I=1,M_DATA
    U_2D_INPUT(I,J) = U_2D_INPUT_2(I,J)*tmp1 +U_2D_INPUT_1(I,J)*tmp2
    V_2D_INPUT(I,J) = V_2D_INPUT_2(I,J)*tmp1 +V_2D_INPUT_1(I,J)*tmp2
    ETA_2D_INPUT(I,J) = ETA_2D_INPUT_2(I,J)*tmp1 +ETA_2D_INPUT_1(I,J)*tmp2
   ENDDO
   ENDDO

120 CONTINUE  ! no more data 

END SUBROUTINE MAP2D_READ_DATA

SUBROUTINE MAP2D_BC
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL
# if defined (PARALLEL)
  USE GLOBAL,ONLY : iista,jjsta   
# endif
                    
  IMPLICIT NONE

! wave components using time-spectra module
    REAL(SP),DIMENSION(NumFreq) :: BB,CC    
    INTEGER :: KK

     ETA_2D_MAP = ZERO
     U_2D_MAP = ZERO
     V_2D_MAP = ZERO


! interpolate data

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       ETA_2D_MAP(I,J) = ETA_2D_MAP(I,J)+(S1(i,j)*ETA_2D_INPUT(nx1(i,j),ny1(i,j)) &
                        +S2(i,j)*ETA_2D_INPUT(nx2(i,j),ny2(i,j)) &
                        +S3(i,j)*ETA_2D_INPUT(nx3(i,j),ny3(i,j))) &
                        /Sc(i,j)

       U_2D_MAP(I,J) =  U_2D_MAP(I,J)+(S1(i,j)*U_2D_INPUT(nx1(i,j),ny1(i,j)) &
                        +S2(i,j)*U_2D_INPUT(nx2(i,j),ny2(i,j)) &
                        +S3(i,j)*U_2D_INPUT(nx3(i,j),ny3(i,j))) &
                        /Sc(i,j)
       V_2D_MAP(I,J) =  V_2D_MAP(I,J)+(S1(i,j)*V_2D_INPUT(nx1(i,j),ny1(i,j)) &
                        +S2(i,j)*V_2D_INPUT(nx2(i,j),ny2(i,j)) &
                        +S3(i,j)*V_2D_INPUT(nx3(i,j),ny3(i,j))) &
                        /Sc(i,j)
     ENDDO
     ENDDO



     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       IF(MASK(I,J)==1)THEN
         ETA(I,J)=ETA_2D_MAP(I,J) +(ETA(I,J)-ETA_2D_MAP(I,J))*SPONGE_2D(I,J)
         U(I,J)=U_2D_MAP(I,J) +(U(I,J)-U_2D_MAP(I,J))*SPONGE_2D(I,J)
         V(I,J)=V_2D_MAP(I,J) +(V(I,J)-V_2D_MAP(I,J))*SPONGE_2D(I,J)
         HU(I,J) = (Depth(I,J)+ETA(I,J))*U(I,J)
         HV(I,J) = (Depth(I,J)+ETA(I,J))*V(I,J)
       ENDIF
     ENDDO
     ENDDO



!print*,'....................'
!print*,V(150,1),V(150,2),V(150,3),V(150,4),V(150,5)


END SUBROUTINE MAP2D_BC

SUBROUTINE MAPPING_COEFF 
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL
# if defined (PARALLEL)
  USE GLOBAL,ONLY : iista,jjsta   
# endif
                    
  IMPLICIT NONE

  REAL(SP) :: area1,area2,area3,x1,x2,x3,y1,y2,y3
  INTEGER :: ii,jj


! 1) calculate x_grid2,y_grid2

! Xco, and Yco


    DO J = 1,Nloc

# if defined (PARALLEL)
	X_grid2(Ibeg,J) = (iista-1)*DX
# else
     X_grid2(Ibeg,J) = ZERO
# endif
     DO I = Ibeg+1,Mloc
       X_grid2(I,J) = X_grid2(I-1,J)+DX
     ENDDO
     DO I = Ibeg-1,Ibeg-Nghost,-1
       X_grid2(I,J) = X_grid2(I+1,J)-DX
     ENDDO

    ENDDO

    DO I = 1,Mloc

# if defined (PARALLEL)
      Y_grid2(I,Jbeg) = (jjsta-1)*DY
# else
     Y_grid2(I,Jbeg) = ZERO
# endif
     DO J = Jbeg+1,Nloc
       Y_grid2(I,J) = Y_grid2(I,J-1)+DY
     ENDDO
     DO J = Jbeg-1,Jbeg-Nghost,-1
       Y_grid2(I,J) = Y_grid2(I,J+1)-DY
     ENDDO
    ENDDO


! 2) s1 s2 s3 ...

        do j=Jbeg,Jend
        do i=Ibeg,Iend
          x1=x_grid2(i,j)
          y1=y_grid2(i,j)
          do jj=1,N_DATA-1
          do ii=1,M_DATA-1

            x2=x_grid1(ii+1,jj)
            y2=y_grid1(ii+1,jj)
            x3=x_grid1(ii,jj+1)
            y3=y_grid1(ii,jj+1)
            area1=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=x_grid1(ii,jj+1)
            y2=y_grid1(ii,jj+1)
            x3=x_grid1(ii,jj)
            y3=y_grid1(ii,jj)
            area2=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=x_grid1(ii,jj)
            y2=y_grid1(ii,jj)
            x3=x_grid1(ii+1,jj)
            y3=y_grid1(ii+1,jj)
            area3=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            if(area1.ge.0.and.area2.ge.0.and.area3.ge.0)then
              nx1(i,j)=ii
              ny1(i,j)=jj
              nx2(i,j)=ii+1
              ny2(i,j)=jj
              nx3(i,j)=ii
              ny3(i,j)=jj+1
              S1(i,j)=area1
              S2(i,j)=area2
              S3(i,j)=area3
               
              x1=x_grid1(nx1(i,j),ny1(i,j))
              y1=y_grid1(nx1(i,j),ny1(i,j))
              x2=x_grid1(nx2(i,j),ny2(i,j))
              y2=y_grid1(nx2(i,j),ny2(i,j))
              x3=x_grid1(nx3(i,j),ny3(i,j))
              y3=y_grid1(nx3(i,j),ny3(i,j))
              Sc(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              goto 110
            endif

            x2=x_grid1(ii+1,jj)
            y2=y_grid1(ii+1,jj)
            x3=x_grid1(ii+1,jj+1)
            y3=y_grid1(ii+1,jj+1)
            area1=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=x_grid1(ii+1,jj+1)
            y2=y_grid1(ii+1,jj+1)
            x3=x_grid1(ii,jj+1)
            y3=y_grid1(ii,jj+1)
            area2=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=x_grid1(ii,jj+1)
            y2=y_grid1(ii,jj+1)
            x3=x_grid1(ii+1,jj)
            y3=y_grid1(ii+1,jj)
            area3=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            if(area1.ge.0.and.area2.ge.0.and.area3.ge.0)then
              nx1(i,j)=ii
              ny1(i,j)=jj+1
              nx2(i,j)=ii+1
              ny2(i,j)=jj
              nx3(i,j)=ii+1
              ny3(i,j)=jj+1
              S1(i,j)=area1
              S2(i,j)=area2
              S3(i,j)=area3
               
              x1=x_grid1(nx1(i,j),ny1(i,j))
              y1=y_grid1(nx1(i,j),ny1(i,j))
              x2=x_grid1(nx2(i,j),ny2(i,j))
              y2=y_grid1(nx2(i,j),ny2(i,j))
              x3=x_grid1(nx3(i,j),ny3(i,j))
              y3=y_grid1(nx3(i,j),ny3(i,j))
              Sc(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              goto 110
            endif

# if defined (PARALLEL)
        IF(MYID==0)  &
         WRITE(*,*) 'Point outside large-model domain, STOP!'
         WRITE(*,*) 'ii,jj : ', II,JJ
        call MPI_FINALIZE ( ier )
# else
         WRITE(*,*) 'Point outside large-model domain, STOP!'   
         WRITE(*,*) 'ii,jj : ', II,JJ    
# endif
         STOP
            
          enddo
          enddo

110       continue

        enddo
        enddo

END SUBROUTINE MAPPING_COEFF 
 
END MODULE ABS_GEN_2D_MODULE

# endif 
! end module
    



