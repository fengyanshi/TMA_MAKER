
  REAL,DIMENSION(:,:),ALLOCATABLE :: Amp_Ser
  REAL,DIMENSION(:),ALLOCATABLE :: Per_Ser,Theta_Ser
  INTEGER :: Ifile,NumFreq,NumDir,I,J
  CHARACTER*80 :: WHAT

  Ifile=1

  OPEN(Ifile,FILE='spectra_file.txt')

! read file
         READ(Ifile,'(A80)')  WHAT ! title

         WRITE(*,*) WHAT
         READ(Ifile,'(A80)')  WHAT ! number of freq  and direction bins
         WRITE(*,*) WHAT

         READ(Ifile,*)  NumFreq,NumDir
           ALLOCATE (Amp_Ser(NumFreq,NumDir),Per_Ser(NumFreq),Theta_Ser(NumDir))

         READ(Ifile,'(A80)')  WHAT ! frequency bins

         WRITE(*,*) WHAT

       DO J=1,NumFreq
          READ(Ifile,*)Per_Ser(J)  ! read in as frequency
!          Per_Ser(J)=1.0/Per_Ser(J)
       ENDDO
         READ(Ifile,'(A80)')  WHAT ! direction bins

         WRITE(*,*) WHAT

       DO I=1,NumDir
          READ(Ifile,*)Theta_Ser(I)
          ! Theta_Ser(I) = Theta_Ser(I)*PI/180.0
       ENDDO
       CLOSE(Ifile)
 
      
       OPEN(Ifile,FILE='spectra_data_1.txt')
         DO I=1,NumDir
           READ(Ifile,*)(Amp_Ser(J,I),J=1,NumFreq)
         ENDDO
       CLOSE(Ifile)


       OPEN(2,FILE='brocch_freq.txt')
       DO J=1,NumFreq
          WRITE(2,*)Per_Ser(J)  
       ENDDO
       CLOSE(2)

       OPEN(2,FILE='brocch_dir.txt')
       DO I=1,NumDir
          WRITE(2,*)Theta_Ser(I)  
       ENDDO
       CLOSE(2)

       OPEN(2,FILE='brocch_data.txt')
       DO I=1,NumDir
          WRITE(2,100)(Amp_Ser(J,I),J=1,NumFreq)
       ENDDO
100    FORMAT(100F12.6)
       CLOSE(2)

       END