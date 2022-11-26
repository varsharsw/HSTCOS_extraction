      SUBROUTINE COSEXTRACT
      CHARACTER LINE*80, ACTPAR*20

      PARAMETER IXBIN = 1
      PARAMETER IYBIN = 1
      PARAMETER IXDIM = 1+16383/IXBIN
      PARAMETER IYDIM = 1+ 1023/IYBIN

      REAL, DIMENSION(IXDIM) :: SPECTRUM, XLAMDA, SPECSMOOTH, QRATIO, RESI, SPEC1, SPEC2,SPECSMOOTH1,SPECSMOOTH2, Q, SP1, SP2
      
      REAL, DIMENSION (IYDIM,IXDIM) :: FIMAGE, RESIDUAL
      REAL, DIMENSION (359,16384) ::PROF
      DIMENSION ROWSUM(1024), Y(1024)
      DIMENSION XR(1024)
      CHARACTER TWOCHAR*3
      REAL DIPOS

      IXMAX = 0
      IYMAX = 0

      FIMAGE = 0
c      CALL coext.py
c      OPEN (UNIT=11,FILE='M33/hst/raw/2dimg/ldy5a3lqq_flt_a_t1.txt', STATUS='OLD')
      OPEN (UNIT=11,FILE='M33/hst/raw/2dimg/ldy5b2vkq_flt_a_t1.txt', STATUS='OLD')
c      OPEN (UNIT=11,FILE='M33/hst/raw/2dimg/laa110gjq_flt_a_t1.txt', STATUS='OLD')
      
      DO IROW=1, 1024
         READ (11, *, err=99) (FIMAGE(IROW,ICOL), ICOL=1,16384)
      ENDDO
      CLOSE (11)
C***  end of reading loop

C***  Extrac spectrum
      SPECTRUM = .0
      DO IROW=410, 460
         DO K = 1, IXDIM
            SPECTRUM(K) = SPECTRUM(K) + FIMAGE(IROW,K)
         ENDDO
      ENDDO

C***  Smooth over plusminus ISMOOTH pixels
      ISMOOTH = 50
      DO K=1, IXDIM
         SPECSMOOTH(K) = 0.
         DO IS = -ISMOOTH, ISMOOTH
            KEXIST = MAX (1, K+IS)        
            KEXIST = MIN (IXDIM,KEXIST) 
            SPECSMOOTH(K) = SPECSMOOTH(K) + SPECTRUM(KEXIST)
         ENDDO
         SPECSMOOTH(K) =SPECSMOOTH(K) / FLOAT(2*ISMOOTH+1)
      ENDDO   
         
      CRPIX1  =    8192.999999999998
      CRVAL1  =    1838.326306296565
      CDELT1  =  0.08001065179671359
      DO K=1, IXDIM
c        XLAMDA(K)=FLOAT(K)
        XLAMDA(K)= (FLOAT(K)-CRPIX1)*CDELT1 + CRVAL1
      ENDDO
      
      KANAL=20
      OPEN (KANAL, FILE='spectrum.plot', STATUS='UNKNOWN')
      CALL PLOTANF (KANAL,'',''
     >        ,'\CENTER\#l# (index)'
     >        ,'\CENTER\Flux'
     >        ,0.,0.,0.,0,0,.0
     >        ,0.,0.,0.,0.,0.,0.
     >        ,XLAMDA,SPECSMOOTH,10000, 5)
      

      OPEN (UNIT=12,FILE='M33/hst/raw/2dimg/prof_tab.txt', STATUS='OLD')
      PROF = 0
      DO IROW=1, 359
         READ (12, *, err=99) (PROF(IROW,ICOL), ICOL=1,16384)
      ENDDO
      CLOSE (12)
C***  end of reading loop

C** re normalization of cross dispersion profile
      DO ICOL=1,16384
          PNORM = .0
          DO IROW=1, 359
             PNORM = PNORM + PROF(IROW,ICOL)
          ENDDO
          DO IROW=1, 359
             PROF(IROW,ICOL) = PROF(IROW,ICOL) / PNORM
          ENDDO
      ENDDO

      DO IROW=1,1024
        XR(IROW)=IROW
        ROWSUM(IROW) = .0
        DO ICOL = 1, 10000
           ROWSUM(IROW) = ROWSUM(IROW) + FIMAGE(IROW,ICOL)
        ENDDO
      ENDDO

      KANAL=2
      OPEN (KANAL, FILE='crossprofile.plot', STATUS='UNKNOWN')
      CALL PLOTANFS (KANAL,'',''
     >        ,''
     >        ,''
     >        ,0.,0.,0.,0,0,.0
     >        ,0.,0.,0.,0.,0.,0.
     >        ,XR,ROWSUM,1024, 'COLOR=2')
      

c     Position of source 1
c      IPOS1=256
c      WRITE (*,*) 'Enter IPOS:'
c      READ (*,*) IPOS1

C***  Offset of source2 from source1
      OFFSET1=245
      OFFSET2=265
      DIPOS=6
      NQ = 111 ! Number of tested qratios
C***  Find optimum QRATIO for each wavelength
      DO IPOS1=OFFSET1,OFFSET2
         DO K = 1, 10000
             DO IQ = 1, NQ
                QRATIO(K) = FLOAT(IQ-1) / FLOAT(NQ-1)
                SPEC1(K) = SPECSMOOTH(K) * QRATIO(K)
                SPEC2(K) = SPECSMOOTH(K) * (1.-QRATIO(K))
                RES=0
                DO IROW=1,359
                   RES = RES + (SPEC1(K) * PROF(IROW,K) + SPEC2(K) * PROF(IROW+DIPOS,K) - FIMAGE(IROW+IPOS1,K))**2
                ENDDO
                IF (IQ .EQ. 1) THEN
                   IQMIN = 1
                   RESMIN = RES
                ELSE
                   IF (RES .LT. RESMIN) THEN
                      IQMIN = IQ
                      RESMIN = RES
                   ENDIF
                ENDIF 
             ENDDO
             QRATIO(K) = FLOAT(IQMIN-1) / FLOAT(NQ-1)
*         PRINT *, DIPOS,Q,RES
         ENDDO
c         CALL PLOTCONS (20, XLAMDA,QRATIO,10000)
         DO K=1, 10000
            SPEC1(K) = SPECTRUM(K) * QRATIO(K) 
            SPEC2(K) = SPECTRUM(K) * (1.-QRATIO(K))
         ENDDO    

C***  Smooth over plusminus ISMOOTH pixels
         DO K=1, IXDIM
            SPECSMOOTH1(K) = 0.
            DO IS = -ISMOOTH, ISMOOTH
               KEXIST = MAX (1, K+IS)        
               KEXIST = MIN (IXDIM,KEXIST) 
               SPECSMOOTH1(K) = SPECSMOOTH1(K) + SPEC1(KEXIST)
            ENDDO
            SPECSMOOTH1(K) =SPECSMOOTH1(K) / FLOAT(2*ISMOOTH+1)
         ENDDO   
C***     Smooth over plusminus ISMOOTH pixels
         DO K=1, IXDIM
            SPECSMOOTH2(K) = 0.
            DO IS = -ISMOOTH, ISMOOTH
               KEXIST = MAX (1, K+IS)        
               KEXIST = MIN (IXDIM,KEXIST) 
               SPECSMOOTH2(K) = SPECSMOOTH2(K) + SPEC2(KEXIST)
            ENDDO
            SPECSMOOTH2(K) =SPECSMOOTH2(K) / FLOAT(2*ISMOOTH+1)
         ENDDO   
         
         WRITE (TWOCHAR, '(I3)') (IPOS1-255)
         
c         CALL PLOTCONS (20, XLAMDA,SPECSMOOTH1,10000,  'COLOR='// TWOCHAR)
c         CALL PLOTCONS (20, XLAMDA,SPECSMOOTH2,10000,  'COLOR='// TWOCHAR)
         
 
c      Q=0.806
         RESIDUAL=FIMAGE
         DO IROW=1,359
            DO K = 1, 10000
               RESIDUAL(IROW+IPOS1,K) = RESIDUAL(IROW+IPOS1,K)-SPEC1(K) *PROF(IROW,K) - SPEC2(K)* PROF(IROW+DIPOS,K)
               ENDDO
         ENDDO
      
C***  Write residual image as big table (for visual inspection as fits file)
c      OPEN (13, FILE='img.dat', STATUS='UNKNOWN')
c      DO IROW=1,1024
c            WRITE (13, '(10000(G14.6,1X))', err=99) (RESIDUAL(IROW,ICOL), ICOL=1,10000)
c      ENDDO 
c      CLOSE(13)
         RESVAL=0
         DO IROW=1,1024
           XR(IROW)=IROW
           ROWSUM(IROW) = .0
           DO ICOL = 1, 10000
              ROWSUM(IROW) = ROWSUM(IROW) + RESIDUAL(IROW,ICOL)
           ENDDO
           RESVAL=RESVAL+ ROWSUM(IROW)**2
         ENDDO

         IF (IPOS1 .EQ. OFFSET1) THEN
              RVALMIN = RESVAL
c              IPOSMIN=IPOS1
c              Q=QRATIO
              SP1=SPECSMOOTH1
              SP2=SPECSMOOTH2
              Y=ROWSUM
         ELSE
              IF (RESVAL .LT. RVALMIN) THEN
                   RVALMIN = RESVAL
c                   IPOSMIN=IPOS1
c                   Q=QRATIO
                   SP1=SPECSMOOTH1
                   SP2=SPECSMOOTH2
                   Y=ROWSUM
c                   PRINT *,IPOS1,MINVAL(ROWSUM)
              ENDIF
         ENDIF 
c         WRITE (*,*) 'TWOCHAR='//TWOCHAR
c         CALL PLOTCONS (2, XR(400),ROWSUM(400),80, 'COLOR='// TWOCHAR)
      ENDDO

      CALL PLOTCONS (20, XLAMDA,SP1,10000, 'COLOR=2')
      CALL PLOTCONS (20, XLAMDA,SP2,10000, 'COLOR=4')
      CALL PLOTCONS (2,XR,Y,1024, 'COLOR=4')
      close(2) 
      CLOSE (20)


      
      STOP 'O.K.'

   99 WRITE (0,*) '***** ERROR when decoding integer number:'
      WRITE (0,*) LINE
      

      END

