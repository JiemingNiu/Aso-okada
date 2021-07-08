      SUBROUTINE NLOKADA_3P(BETA, LAMBDA, MU,
     *  NSTA, NCHN, NY, 
     *  SXX, SYY, SZZ, 
     *  RXX, RYY, RZZ1, RZZ2, Y_O, DY_O, 
     *  X_M, PHI_M, THETA_M, RMS_M)
      IMPLICIT NONE

      INTEGER          NSTA, NCHN, NX, NY, NW, NYR
      DOUBLE PRECISION LAMBDA, MU, BETA, SXX, SYY, SZZ
      DOUBLE PRECISION PHI_M, THETA_M, RMS_M
      DOUBLE PRECISION RXX(NSTA), RYY(NSTA)
      DOUBLE PRECISION RZZ1(NSTA), RZZ2(NSTA)
      DOUBLE PRECISION Y_O(NY), DY_O(NY)

      INTEGER          NTHETA, NPHI
      INTEGER          IX, IY, IPHI, ITHETA, ISTA, J
      DOUBLE PRECISION ASTEP, DEGTORAD, PHI, RMS_M_, S0, V 
      PARAMETER        (NX=4, NW=20)
      PARAMETER        (ASTEP=1.0, DEGTORAD=1.745329252E-02)

      DOUBLE PRECISION Y_ON(NY+NX), DY_ON(NY), Y_MN(NY+NX)
      DOUBLE PRECISION X_M(NX), X_M_(NX)
      DOUBLE PRECISION A_M(NY+NX,NX), A_M_(NY+NX, NX)

!FROM OKADA'S SUBROUTINE DC3D0:
      INTEGER          IRET, INFO
      REAL*4           ALPHA,POT1,POT2,POT3,POT4,X,Y,Z,DEPTH,DIP
      REAL*4           UXS,UYS,UZS,UXXS,UYXS,UZXS,
     &  UXYS,UYYS,UZYS,UXZS,UYZS,UZZS 
      REAL*4           UXD,UYD,UZD,UXXD,UYXD,UZXD,
     &  UXYD,UYYD,UZYD,UXZD,UYZD,UZZD 
      REAL*4           UXT,UYT,UZT,UXXT,UYXT,UZXT,
     &  UXYT,UYYT,UZYT,UXZT,UYZT,UZZT 
      REAL*4           UXE,UYE,UZE,UXXE,UYXE,UZXE,
     &  UXYE,UYYE,UZYE,UXZE,UYZE,UZZE 
      DOUBLE PRECISION CSST, SSST, CS2ST, SS2ST
      DOUBLE PRECISION RXX_M(NSTA), RYY_M(NSTA), WORK(NW)
      
      !PRINT *, 'SXX=',SXX,';SYY=',SYY,';SZZ=',SZZ
      DO IY=1,NY
        DY_ON(IY) = DY_O(IY) / ABS(Y_O(IY))
        !WRITE(*,*) 'Y_O(IY)=', IY, Y_O(IY)
      END DO
      !
      S0 = 10.0 * DY_ON(1)
      DO IY=1,NY
        DY_ON(IY) = DY_ON(IY) / S0 
        !WRITE(*,*) 'DY_ON(IY)=', IY, DY_ON(IY)
      END DO
      !
      DO IY=1,NY
        Y_ON(IY) = Y_O(IY) / ABS(Y_O(IY)) / DY_ON(IY)
      END DO
      DO IY=1,NX
        Y_ON(IY+NY) = 0.0
      END DO

      ALPHA = SNGL((LAMBDA+MU)/(LAMBDA+2.0*MU))
      POT1 = SNGL(1.0E15/MU)
      POT2 = SNGL(1.0E15/MU)
      POT3 = SNGL(1.0E15/LAMBDA)
      POT4 = SNGL(1.0E15/MU)
      !WRITE(*,*), 'POTS=', POT1, POT2, POT3, POT4
 
      NTHETA = int(90.0/ASTEP)
      NPHI = int(360.0/ASTEP)
      !NTHETA = 1
      !NPHI = 1

      RMS_M = 9999.9
      DO IPHI=1,NPHI
        PHI = ASTEP * (IPHI-1) * DEGTORAD
        CSST=DCOS(PHI)
        SSST=DSIN(PHI)

        DO ISTA=1,NSTA
          RXX_M(ISTA)=SNGL((RXX(ISTA)-SXX)*SSST+(RYY(ISTA)-SYY)*CSST)
          RYY_M(ISTA)=SNGL(-(RXX(ISTA)-SXX)*CSST+(RYY(ISTA)-SYY)*SSST)
        END DO
        
        DO ITHETA=1,NTHETA
          DIP = SNGL(ASTEP * (ITHETA-1))
          DO ISTA=1,NSTA
            J = (ISTA-1)*3
            X=SNGL(RXX_M(ISTA))
            Y=SNGL(RYY_M(ISTA))
            !Z=SNGL(RZZ1(ISTA)-RZZ2(ISTA))
            !WRITE(*,*) '-->', ISTA, X, Y
! calculate displacement at the free surface: Z=0.0
            DEPTH=SNGL(RZZ1(ISTA)+SZZ)
            CALL DC3D0(ALPHA,X,Y,0.0,DEPTH,DIP,POT1,0.0,0.0,0.0,
     *       UXS,UYS,UZS,UXXS,UYXS,UZXS,UXYS,UYYS,
     *       UZYS,UXZS,UYZS,UZZS,IRET)
            CALL DC3D0(ALPHA,X,Y,0.0,DEPTH,DIP,0.0,POT2,0.0,0.0,
     *       UXD,UYD,UZD,UXXD,UYXD,UZXD,UXYD,UYYD,
     *       UZYD,UXZD,UYZD,UZZD,IRET)
            CALL DC3D0(ALPHA,X,Y,0.0,DEPTH,DIP,0.0,0.0,POT3,0.0,
     *       UXT,UYT,UZT,UXXT,UYXT,UZXT,UXYT,UYYT,
     *       UZYT,UXZT,UYZT,UZZT,IRET)
            CALL DC3D0(ALPHA,X,Y,0.0,DEPTH,DIP,0.0,0.0,0.0,POT4,
     *       UXE,UYE,UZE,UXXE,UYXE,UZXE,UXYE,UYYE,
     *       UZYE,UXZE,UYZE,UZZE,IRET)
            A_M(J+1, 1) = DBLE(UZS)
            A_M(J+1, 2) = DBLE(UZD)
            A_M(J+1, 3) = DBLE(UZT)
            A_M(J+1, 4) = DBLE(UZE)

! calculate tilt at the borehole depth: Z<0.0
            Z = -SNGL(RZZ2(ISTA))
            CALL DC3D0(ALPHA,X,Y,Z,DEPTH,DIP,POT1,0.0,0.0,0.0,
     *       UXS,UYS,UZS,UXXS,UYXS,UZXS,UXYS,UYYS,
     *       UZYS,UXZS,UYZS,UZZS,IRET)
            CALL DC3D0(ALPHA,X,Y,Z,DEPTH,DIP,0.0,POT2,0.0,0.0,
     *       UXD,UYD,UZD,UXXD,UYXD,UZXD,UXYD,UYYD,
     *       UZYD,UXZD,UYZD,UZZD,IRET)
            CALL DC3D0(ALPHA,X,Y,Z,DEPTH,DIP,0.0,0.0,POT3,0.0,
     *       UXT,UYT,UZT,UXXT,UYXT,UZXT,UXYT,UYYT,
     *       UZYT,UXZT,UYZT,UZZT,IRET)
            CALL DC3D0(ALPHA,X,Y,Z,DEPTH,DIP,0.0,0.0,0.0,POT4,
     *       UXE,UYE,UZE,UXXE,UYXE,UZXE,UXYE,UYYE,
     *       UZYE,UXZE,UYZE,UZZE,IRET)
 
            A_M(J+2, 1) = -(DBLE(UZXS)*SSST-DBLE(UZYS)*CSST)
            A_M(J+2, 2) = -(DBLE(UZXD)*SSST-DBLE(UZYD)*CSST)
            A_M(J+2, 3) = -(DBLE(UZXT)*SSST-DBLE(UZYT)*CSST)
            A_M(J+2, 4) = -(DBLE(UZXE)*SSST-DBLE(UZYE)*CSST)

            A_M(J+3, 1) = -(DBLE(UZXS)*CSST+DBLE(UZYS)*SSST)
            A_M(J+3, 2) = -(DBLE(UZXD)*CSST+DBLE(UZYD)*SSST)
            A_M(J+3, 3) = -(DBLE(UZXT)*CSST+DBLE(UZYT)*SSST)
            A_M(J+3, 4) = -(DBLE(UZXE)*CSST+DBLE(UZYE)*SSST)

          END DO
          !			CHOLESKY FACTORISATION
          !WRITE(*,*) 'PHI=',phi,'THETA=', DIP
          !WRITE(*,*) 'Before normliasing:'
          !WRITE(*,*) 'A_M='
          !DO IY=1,NY
          !  PRINT *, A_M(IY,:)
          !END DO
          DO IY=1,NY
            DO IX=1,NX
              A_M(IY, IX) = A_M(IY, IX) / ABS(Y_O(IY)) / DY_ON(IY)
            END DO
            Y_MN(IY) = Y_ON(IY)
          END DO
          !
          DO IY=1,NX
            Y_MN(IY+NY) = Y_ON(IY+NY)
            DO IX=1,NX
              IF (IX.EQ.IY) THEN
                A_M(IY+NY, IX) = BETA
              ELSE
                A_M(IY+NY, IX) = 0.0
              END IF
            END DO
          END DO
          !
          DO IY=1,NY+NX
            DO IX=1,NX
              A_M_(IY, IX) = A_M(IY, IX)
            END DO
          END DO
          !
          !WRITE(*,*) 'After normliasing:'
          !WRITE(*,*) 'A_M='
          !DO IY=1,NY+NX
          !  PRINT *, A_M(IY,:)
          !  WRITE(*,*) 'Y_ON(IY)=', IY, Y_MN(IY)
          !END DO
          !LEAST SQUARE
          ! dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
          CALL DGELS('N', NY+NX, NX, 1, A_M_, NY+NX, Y_MN, NY+NX, 
     *     WORK, NW, INFO)
          DO IX=1,NX
            X_M_(IX) = Y_MN(IX)
            !WRITE(*,*) 'X_M=', IX, X_M(IX)
          END DO

          RMS_M_ = 0.0
          DO IY=1,NY
            V = 0.0
            DO IX=1,NX
              V = V + A_M(IY, IX) * X_M_(IX)
            END DO
            V = V - Y_ON(IY)
            !PRINT *, IY, V, Y_ON(IY)
            RMS_M_ = RMS_M_ + V * V 
          END DO
          
          IF(RMS_M_.LE.RMS_M) THEN
            RMS_M = RMS_M_
            PHI_M = PHI / DEGTORAD
            THETA_M = DBLE(DIP)
            DO IX=1,NX
              X_M(IX) = X_M_(IX)
            END DO
          END IF
        END DO
      END DO
      END SUBROUTINE NLOKADA_3P 

