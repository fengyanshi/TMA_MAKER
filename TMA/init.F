!------------------------------------------------------------------------------------
!
!      FILE init.F
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
!    ALLOCATE_VARIABLES is a subroutine to allocate variables
!
!    HISTORY:
!    05/01/2010 Fengyan Shi 
!    09/26/2013 Babak Tehranirad, added spacially varying cd
!    08/06/2015 Choi added V4xL, U4yL, V4xL, U4yL, WAVEMAKER_VIS
!
!-------------------------------------------------------------------------------------
SUBROUTINE ALLOCATE_VARIABLES
     USE GLOBAL

! coordinate for cartesian only
     ALLOCATE (Xco(Mloc),Yco(Nloc))
! allocate variables
     ALLOCATE(DelxU(Mloc,Nloc),DelxHU(Mloc,Nloc),DelxV(Mloc,Nloc),DelxEtar(Mloc,Nloc),&
              DelyU(Mloc,Nloc),DelyHV(Mloc,Nloc),DelyV(Mloc,Nloc),DelyEtar(Mloc,Nloc),&
              DelxHV(Mloc,Nloc),DelyHU(Mloc,Nloc), &
! U V HU H in x-direction
              UxL(Mloc1,Nloc),UxR(Mloc1,Nloc),VxL(Mloc1,Nloc),VxR(Mloc1,Nloc),&
              HUxL(Mloc1,Nloc),HUxR(Mloc1,Nloc),HVxL(Mloc1,Nloc),HVxR(Mloc1,Nloc), &
              HxL(Mloc1,Nloc),HxR(Mloc1,Nloc), &
! U V HV H in y-direction
              UyL(Mloc,Nloc1),UyR(Mloc,Nloc1),VyL(Mloc,Nloc1),VyR(Mloc,Nloc1),&
              HVyL(Mloc,Nloc1),HVyR(Mloc,Nloc1),HUyL(Mloc,Nloc1),HUyR(Mloc,Nloc1), &
              HyL(Mloc,Nloc1),HyR(Mloc,Nloc1), &
! cross-derivatives
              Uxy(Mloc,Nloc),Vxy(Mloc,Nloc),DUxy(Mloc,Nloc),DVxy(Mloc,Nloc), &
! second-derivatives
              Uxx(Mloc,Nloc),Vyy(Mloc,Nloc),DUxx(Mloc,Nloc),DVyy(Mloc,Nloc), &
! 1st-derivatives
              Ux(Mloc,Nloc),Uy(Mloc,Nloc),Vx(Mloc,Nloc),Vy(Mloc,Nloc), &
              DUx(Mloc,Nloc),DUy(Mloc,Nloc),DVx(Mloc,Nloc),DVy(Mloc,Nloc), &
              ETAT(Mloc,Nloc),ETAx(Mloc,Nloc),ETAy(Mloc,Nloc), &
              ETATx(Mloc,Nloc),ETATy(Mloc,Nloc), &
! time-derivatives
              U0(Mloc,Nloc),V0(Mloc,Nloc),Ut(Mloc,Nloc),Vt(Mloc,Nloc),&
              Utx(Mloc,Nloc),Vty(Mloc,Nloc),Utxx(Mloc,Nloc),Utxy(Mloc,Nloc),&
              Vtxy(Mloc,Nloc),Vtyy(Mloc,Nloc),&
              DUtxx(Mloc,Nloc),DUtxy(Mloc,Nloc),&
              DVtxy(Mloc,Nloc),DVtyy(Mloc,Nloc),DUtx(Mloc,Nloc),DVty(Mloc,Nloc),&
! P Q Eta, Fx, Fy
              PL(Mloc1,Nloc),PR(Mloc1,Nloc),QL(Mloc,Nloc1),QR(Mloc,Nloc1), &
              FxL(Mloc1,Nloc),FxR(Mloc1,Nloc),FyL(Mloc,Nloc1),FyR(Mloc,Nloc1), &
              GxL(Mloc1,Nloc),GxR(Mloc1,Nloc),GyL(Mloc,Nloc1),GyR(Mloc,Nloc1), &
              EtaRxL(Mloc1,Nloc),EtaRxR(Mloc1,Nloc), &
              EtaRyL(Mloc,Nloc1),EtaRyR(Mloc,Nloc1), &
! sponge
              SPONGE(Mloc,Nloc), SpongeMaker(Mloc,Nloc), &
! original variables at notes
              Fx(Mloc1,Nloc),Fy(Mloc,Nloc1),&
              U(Mloc,Nloc),V(Mloc,Nloc), HU(Mloc,Nloc),HV(Mloc,Nloc),&
              Gx(Mloc1,Nloc),Gy(Mloc,Nloc1), &
              P(Mloc1,Nloc),Q(Mloc,Nloc1), &
              SxL(Mloc1,Nloc),SxR(Mloc1,Nloc), &
              SyL(Mloc,Nloc1),SyR(Mloc,Nloc1),SourceX(Mloc,Nloc), &
              SourceY(Mloc,Nloc), &
! others
              Umean(Mloc,Nloc),Vmean(Mloc,Nloc),ETAmean(Mloc,Nloc),&
              Usum(Mloc,Nloc),Vsum(Mloc,Nloc),ETAsum(Mloc,Nloc), &
              nu_smg(Mloc,Nloc), &
              Num_Zero_Up(Mloc,Nloc), &
              WaveHeightRMS(Mloc,Nloc),  &
              WaveHeightAve(Mloc,Nloc),  &
              Emax(Mloc,Nloc),  &
              Emin(Mloc,Nloc), &
              HrmsSum(Mloc,Nloc), &
              HavgSum(Mloc,Nloc), &
	        !ykchoi
	        ETA2sum(Mloc,Nloc), ETA2mean(Mloc,Nloc), &
			SigWaveHeight(Mloc,Nloc),  &
# if defined (CARTESIAN)
              U4xL(Mloc1,Nloc),U4xR(Mloc1,Nloc),&
              V4yL(Mloc,Nloc1),V4yR(Mloc,Nloc1), &
	! ykchoi added V4xL and U4yL (08/06/15)
			V4xL(Mloc1,Nloc),V4xR(Mloc1,Nloc),&  
			U4yL(Mloc,Nloc1),U4yR(Mloc,Nloc1) & 
# else
# if defined (ZALPHA)
              U4xL(Mloc1,Nloc),U4xR(Mloc1,Nloc),&
              V4yL(Mloc,Nloc1),V4yR(Mloc,Nloc1), &
        ! ykchoi added V4xL and U4yL (08/06/15)
	      V4xL(Mloc1,Nloc),V4xR(Mloc1,Nloc),&  
	      U4yL(Mloc,Nloc1),U4yR(Mloc,Nloc1), &  
# endif
              Dx(Mloc,Nloc),Dy(Mloc,Nloc), Coriolis(Mloc,Nloc), &
              SlopeX(Mloc,Nloc),SlopeY(Mloc,Nloc), &
              Lat_theta(Mloc,Nloc) &
# endif
              )
      ALLOCATE(Depth(Mloc,Nloc),H(Mloc,Nloc),&
               Depthx(Mloc1,Nloc),Depthy(Mloc,Nloc1), &
               MASK(Mloc,Nloc),DepthNode(Mloc1,Nloc1), &
               MASK_STRUC(Mloc,Nloc),MASK9(Mloc,Nloc), &
               tmp4preview(Mloc,Nloc),Int2Flo(Mloc,Nloc),&
               Cd(Mloc,Nloc),CD_breakwater(Mloc,Nloc) &
              )
! updating variables
      ALLOCATE(Eta(Mloc,Nloc),Eta0(Mloc,Nloc), &
               Ubar0(Mloc,Nloc),Vbar0(Mloc,Nloc),&
               Ubar(Mloc,Nloc),Vbar(Mloc,Nloc))
! dispersion updating variables
# if defined (CARTESIAN)
      ALLOCATE(U4(Mloc,Nloc),V4(Mloc,Nloc),U1p(Mloc,Nloc), & 
               V1p(Mloc,Nloc),U1pp(Mloc,Nloc),V1pp(Mloc,Nloc),&
               U2(Mloc,Nloc),V2(Mloc,Nloc),U3(Mloc,Nloc),V3(Mloc,Nloc))
# else
# if defined (ZALPHA)
      ALLOCATE(U4(Mloc,Nloc),V4(Mloc,Nloc) & 
               )
# endif
      ALLOCATE(U1p(Mloc,Nloc), & 
               V1p(Mloc,Nloc) )
# endif
# if defined (ITERATION)
      ALLOCATE(Ui(Mloc,Nloc), Vi(Mloc,Nloc),ETAi(Mloc,Nloc),&
        UbarOld(Mloc,Nloc),VbarOld(Mloc,Nloc),EtaOld(Mloc,Nloc) )
# endif

  ! HeightMax will be used not only in output but also meteo module
        ALLOCATE(HeightMax(Mloc,Nloc))
        HeightMax=ZERO

      ALLOCATE(WaveMaker_Mass(Mloc,Nloc))
      WaveMaker_Mass = ZERO

      IF(WAVEMAKER(1:7)=='WK_TIME')THEN
        ALLOCATE(WAVE_COMP(NumWaveComp,3),Beta_genS(NumWaveComp),D_genS(NumWaveComp) )
      ENDIF

      IF(VISCOSITY_BREAKING.OR.SHOW_BREAKING)THEN
       ALLOCATE(AGE_BREAKING(Mloc,Nloc))
       ALLOCATE(nu_break(Mloc,Nloc))
       nu_break=nu_bkg
      ENDIF
	
      IF(WAVEMAKER_VIS)THEN
       ALLOCATE(nu_break(Mloc,Nloc))
       nu_break=ZERO
      ENDIF

      IF(DIFFUSION_SPONGE)THEN
       ALLOCATE(nu_sponge(Mloc,Nloc))
       nu_sponge=ZERO
      ENDIF

      IF(OUT_Hmin)THEN
        ALLOCATE(HeightMin(Mloc,Nloc))
        HeightMin=ZERO
      ENDIF
      IF(OUT_Umax)THEN
        ALLOCATE(VelocityMax(Mloc,Nloc))
        VelocityMax=ZERO
      ENDIF
      IF(OUT_VORmax)THEN
        ALLOCATE(VorticityMax(Mloc,Nloc))
        VorticityMax=ZERO
      ENDIF
      IF(OUT_MFmax)THEN
        ALLOCATE(MomentumFluxMax(Mloc,Nloc))
        MomentumFluxMax=ZERO
      ENDIF
      

END SUBROUTINE ALLOCATE_VARIABLES

!-------------------------------------------------------------------------------------
!
!    INITIALIZATION is subroutine for initialization
!
!    HISTORY:
!    05/01/2010 Fengyan Shi
!    09/26/2013 Babak Tehranirad, added varying cd
!
!-------------------------------------------------------------------------------------
SUBROUTINE INITIALIZATION
     USE GLOBAL
     USE INPUT_READ
# if defined(INTEL)
     USE IFPORT
# endif
     IMPLICIT NONE
     INTEGER :: VTYPE
     CHARACTER(LEN=80) :: WHAT
# if defined (PARALLEL)
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: VarGlob
     REAL(SP) :: myvar_tmp
# endif

# if !defined (CARTESIAN)
     REAL(SP)::Dphi_r,Dtheta_r

! for spherical coordinates and Coriolis parameter

   IF(StretchGrid) THEN
# if defined (PARALLEL)
       CALL GetFile(DX_FILE,Dx)
       CALL GetFile(DY_FILE,Dy)
       CALL GetFile(Coriolis_FILE,Coriolis)
# else
      OPEN(1,FILE=TRIM(DX_FILE))
       DO J=Jbeg,Jend
         READ(1,*)(Dx(I,J),I=Ibeg,Iend)
       ENDDO
      CLOSE(1)

      OPEN(1,FILE=TRIM(DY_FILE))
       DO J=Jbeg,Jend
         READ(1,*)(Dy(I,J),I=Ibeg,Iend)
       ENDDO
      CLOSE(1)

      OPEN(1,FILE=TRIM(Coriolis_FILE))
       DO J=Jbeg,Jend
         READ(1,*)(Coriolis(I,J),I=Ibeg,Iend)
       ENDDO
      CLOSE(1)

! ghost cells
     DO I=Ibeg,Iend
       DO J=1,Nghost
        Dx(I,J)=Dx(I,Jbeg)
        Dy(I,J)=Dy(I,Jbeg)
        Coriolis(I,J)=Coriolis(I,Jbeg)
       ENDDO
       DO J=Jend+1,Nloc
        Dx(I,J)=Dx(I,Jend)
        Dy(I,J)=Dy(I,Jend)
        Coriolis(I,J)=Coriolis(I,Jend)
       ENDDO
     ENDDO

     DO J=1,Nloc
       DO I=1,Nghost
        Dx(I,J)=Dx(Ibeg,J)
        Dy(I,J)=Dy(Ibeg,J)
        Coriolis(I,J)=Coriolis(Ibeg,J)
       ENDDO
       DO I=Iend+1,Mloc
        Dx(I,J)=Dx(Iend,J)
        Dy(I,J)=Dy(Iend,J)
        Coriolis(I,J)=Coriolis(Iend,J)
       ENDDO
     ENDDO  

# endif

    ELSE ! for stretch grid

     Dphi_r=Dphi*pi/180.0_SP
     Dtheta_r=Dtheta*pi/180.0_SP

     Do J=1,Nloc
     Do I=1,Mloc
# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!       Lat_theta(I,J)=Lat_South*pi/180.0_SP-Nghost*Dtheta_r &
!                       +(npy*Nglob/py+J-1)*Dtheta_r
       Lat_theta(I,J)=Lat_South*pi/180.0_SP-Nghost*Dtheta_r &
                       +((jjsta - 1)+J-1)*Dtheta_r
!---ykchoi Jan/23/2018]
# else
       Lat_theta(I,J)=Lat_South*pi/180.0_SP-Nghost*Dtheta_r &
                       +(J-1)*Dtheta_r
# endif
       Dx(I,J) = R_earth*Dphi_r*COS(Lat_theta(I,J))
       Dy(I,J) = R_earth*Dtheta_r
       Coriolis(I,J) = pi*SIN(Lat_theta(I,J)) / 21600.0_SP
     ENDDO
     ENDDO
    ENDIF ! stretchgrid

# endif

! parameter kappa for order of MUSCL
     IF(HIGH_ORDER(1:3)=='SEC')THEN
      Kappa = -1.0_SP
     ELSE
      Kappa = 1.0_SP/3.0_SP
     ENDIF

! set zeros
     T_sum = ZERO
     Umean = ZERO
     Vmean = ZERO
     ETAmean = ZERO
     nu_smg = ZERO
     Num_Zero_Up = 0
     WaveHeightRMS = ZERO
     WaveHeightAve =ZERO 
     Emax = ZERO
     Emin = ZERO
     HrmsSum = ZERO
     HavgSum = ZERO
     DelxU=0.0_SP
     DelxHU=0.0_SP
     DelxV=0.0_SP
     DelxEtar=0.0_SP
     DelyU=0.0_SP
     DelyHV=0.0_SP
     DelyV=0.0_SP
     DelyEtar=0.0_SP
     DelxHV=0.0_SP
     DelyHU=0.0_SP
     UxL=0.0_SP
     UxR=0.0_SP
     VxL=0.0_SP
     VxR=0.0_SP
     HUxL=0.0_SP
     HUxR=0.0_SP
     HVxL=0.0_SP
     HVxR=0.0_SP
     HxL=0.0_SP
     HxR=0.0_SP
     UyL=0.0_SP
     UyR=0.0_SP
     VyL=0.0_SP
     VyR=0.0_SP
     HVyL=0.0_SP
     HVyR=0.0_SP
     HUyL=0.0_SP
     HUyR=0.0_SP
     HyL=0.0_SP
     HyR=0.0_SP
# if defined (CARTESIAN)
     U4xL=ZERO
     U4xR=ZERO
     V4yL=ZERO
     V4yR=ZERO
# endif
     Uxy=ZERO
     Vxy=ZERO
     DUxy=ZERO
     DVxy=ZERO
     Uxx=ZERO
     Vyy=ZERO
     DUxx=ZERO
     DVyy=ZERO 
     U0=ZERO
     V0=ZERO
     Ut=ZERO
     Vt=ZERO
     Utx=ZERO
     Vty=ZERO
     Utxx=ZERO
     Utxy=ZERO
     Vtxy=ZERO
     Vtyy=ZERO
     DUtxx=ZERO
     DUtxy=ZERO
     DVtxy=ZERO
     DVtyy=ZERO
     DUtx=ZERO
     DVty=ZERO    
     PL=0.0_SP
     PR=0.0_SP
     QL=0.0_SP
     QR=0.0_SP
     FxL=0.0_SP
     FxR=0.0_SP
     FyL=0.0_SP
     FyR=0.0_SP
     GxL=0.0_SP
     GxR=0.0_SP
     GyL=0.0_SP
     GyR=0.0_SP
     SxL=0.0_SP
     SxR=0.0_SP
     SyL=0.0_SP
     SyR=0.0_SP
! original variables
     Ubar=0.0_SP
     Vbar=0.0_SP
     Ubar0=0.0_SP
     Vbar0=0.0_SP
     U=0.0_SP
     V=0.0_SP
     HU=0.0_SP
     HV=0.0_SP
     Fx=0.0_SP
     Fy=0.0_SP
     Gx=0.0_SP
     Gy=0.0_SP
     P=0.0_SP
     Q=0.0_SP
     U1p=ZERO
     V1p=ZERO
# if defined (CARTESIAN)
     U4=ZERO
     V4=ZERO
     U1pp=ZERO
     V1pp=ZERO
     U2=ZERO
     V2=ZERO
     U3=ZERO
     V3=ZERO
# endif
     Depth=10.0_SP
     DepthNode=10.0_SP
     H=0.0_SP
     Eta=0.0_SP
     SourceX=0.0_SP
     SourceY=0.0_SP
     PLOT_COUNT=0.0_SP
     PLOT_COUNT_STATION=0.0_SP
     HOTSTART_COUNT=ZERO
     MASK=1
     MASK_STRUC=1
     SCREEN_COUNT=ZERO
     SPONGE=1.0_SP
     SpongeMaker=1.0_SP

! coordinate for cartesian only
! Xco, and Yco


     CALL WAVEMAKER_INITIALIZATION
   
   


END SUBROUTINE INITIALIZATION


!-------------------------------------------------------------------------------------
!
!    INITIAL_UVZ is subroutine of given initial u v and eta 
!
!    HISTORY:
!    02/03/2011 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE INITIAL_UVZ
      USE GLOBAL
      IMPLICIT NONE
      REAL(SP),DIMENSION(Mloc,Nloc) :: rMASK
      LOGICAL :: FILE_EXIST

IF(.NOT.NO_UV_FILE)THEN

  INQUIRE(FILE=TRIM(U_FILE),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0) WRITE(*,*) TRIM(U_FILE), ' specified in input.txt but does not exist. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(U_FILE), ' specified in input.txt but does not exist. STOP'
    STOP
# endif
  ENDIF


  INQUIRE(FILE=TRIM(V_FILE),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0) WRITE(*,*) TRIM(V_FILE), ' specified in input.txt but does not exist. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(V_FILE), ' specified in input.txt but does not exist. STOP'
    STOP
# endif
  ENDIF

ENDIF

  INQUIRE(FILE=TRIM(ETA_FILE),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0) WRITE(*,*) TRIM(ETA_FILE), ' specified in input.txt but does not exist. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(ETA_FILE), ' specified in input.txt but does not exist. STOP'
    STOP
# endif
  ENDIF

# if defined (PARALLEL)
IF(.NOT.NO_UV_FILE)THEN
      CALL GetFile(U_FILE,U)
      CALL GetFile(V_FILE,V)
ELSE
      U=ZERO
      V=ZERO
ENDIF

      CALL GetFile(ETA_FILE,ETA)

      IF(.NOT.NO_MASK_FILE)THEN
      CALL GetFile(MASK_FILE,rMASK)
      ENDIF

# else

IF(.NOT.NO_UV_FILE)THEN
      OPEN(1,FILE=TRIM(U_FILE))
       DO J=Nghost+1,Nloc-Nghost
        READ(1,*)(U(I,J),I=Nghost+1,Mloc-Nghost)
       ENDDO
      CLOSE(1)

      OPEN(1,FILE=TRIM(V_FILE))
       DO J=Nghost+1,Nloc-Nghost
        READ(1,*)(V(I,J),I=Nghost+1,Mloc-Nghost)
       ENDDO
      CLOSE(1)
ELSE
      U=ZERO
      V=ZERO
ENDIF

      OPEN(1,FILE=TRIM(ETA_FILE))
       DO J=Nghost+1,Nloc-Nghost
        READ(1,*)(ETA(I,J),I=Nghost+1,Mloc-Nghost)
       ENDDO
      CLOSE(1)

     IF(.NOT.NO_MASK_FILE)THEN
      OPEN(1,FILE=TRIM(MASK_FILE))
       DO J=Nghost+1,Nloc-Nghost
        READ(1,*)(rMASK(I,J),I=Nghost+1,Mloc-Nghost)
       ENDDO
      CLOSE(1)
     ENDIF

# endif

     IF(.NOT.NO_MASK_FILE)THEN
      MASK=INT(rMASK)
     ELSE
      MASK = 1
     ENDIF

END SUBROUTINE INITIAL_UVZ




