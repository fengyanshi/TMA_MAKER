!------------------------------------------------------------------------------------
!
!      FILE main.F
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
Program FUNWAVE_TVD
!-------------------------------------------------------------------------------------
!    VERSION 3.0
!
!    MAIN - READ_INPUT
!         - INDEX
!         - ALLOCATE_VARIABLES
!         - INITIALIZATION 
!         DO LOOP
!             - VARIABLE_UPDATING
!             - EXCHANGE
!             - ESTIMATE_DT
!             - RUNGE-KUTTA
!             - DISPERSION
!             - FLUXES          
!             - SourceTerms      
!             - ESTIMATE_HUV
!             - EXCHANGE 
!             - FLUXES                 
!             - SourceTerms again
!
!             - STATISTICS
!             - PREVIEW
!          ENDDO LOOP
!-------------------------------------------------------------------------------------
! ** OPEN FILES **
!  (1): read input, (2): output, (3): log, (4): !write/read hotstart
!-------------------------------------------------------------------------------------
! ** HOT START DATA **
!   NOTE: read input.txt first, if HOT_START, then read  
        ! -- dimension
! Mloc,Nloc,Mloc1,Nloc1
! Nghost
! Ibeg,Iend,Jbeg,Jend,Iend1,Jend1
!   NOTE: need to confirm if the saved data is consistent with input.txt
        ! -- time
! TIME
! TOTAL_TIME
! PLOT_INTV
! PLOT_COUNT
! SCREEN_INTV
! SCREEN_COUNT
! HOTSTART_INTV
! ICOUNT
        ! spacing
! DX,DY
        ! -- physics
! DISPERSION
! Gamma1
! a1,a2,b1,b2
! SWE_ETA_DEP
        ! -- numerics
! Time_Scheme
! HIGH_ORDER
! CONSTR
! CFL
! FroudeCap
! DISP_TIME_LEFT
        ! -- wet-dry
! MinDepth,MinDepthfrc

        ! -- depth
! DEPTH
! DEPTHx
! DEPTHy
        ! variables
! U
! V
! if (.NOT.DISP_TIME_LEFT)THEN
! U0
! V0
! endif
! Ubar
! Vbar
! ETA 
! H
! MASK
! MASK9
! MAST_STRUC
!
       ! -- wavemaker
! if (WAVEMAKER is WK_IRR)
! turns out the data for Cm Sm too large, calculate it when hotstart
!
! if (WAVEMAKER is WK_REG)
! D_gen
! Beta_gen
! rlamda
! 
!
!-------------------------------------------------------------------------------------
     USE GLOBAL
# if defined (VESSEL)
     USE VESSEL_MODULE
# endif
# if defined (SEDIMENT)
     USE SEDIMENT_MODULE
# endif

# if defined (METEO)
     USE METEO_MODULE
# endif
# if defined (TRACKING)
     USE TRACER
# endif

     IMPLICIT NONE

!     INTEGER::ISTAGE ! moved to mod_global 09/12/2017
!     REAL(SP) :: tbegin,tend  ! moved to mod_global 07/29/2016

# if defined (PARALLEL)
     CALL MPI_INIT ( ier )
# endif

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'MPI_INIT'     
#    endif
# endif

     CALL READ_INPUT

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'READ_INPUT'     
#    else   
       WRITE(*,*) 'READ_INPUT'
#    endif
# endif

	![ykchoi(14.12.24.)
	IF(INI_UVZ)THEN
        TIME=HotStartTime
      ENDIF
	!ykchoi(14.12.24.)]


# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'INDEX'     
#    else   
       WRITE(*,*) 'INDEX'
#    endif
# endif

! allocate variables
     CALL ALLOCATE_VARIABLES

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'ALLOCATION'     
#    else   
       WRITE(*,*) 'ALLOCATION'
#    endif
# endif

     CALL INITIALIZATION
      
# if defined (PARALLEL)
     call MPI_FINALIZE ( ier )
# endif

END PROGRAM FUNWAVE_TVD

