      PROGRAM NLOPT_MAIN
      IMPLICIT NONE
      INCLUDE          "mpif.h"
      INTEGER          MAX_STA, NCHN, NXSTEP, NZSTEP
      DOUBLE PRECISION XSTART, ZSTART, XSTEP, ZSTEP
      DOUBLE PRECISION VP, VS, RHO
      !PARAMETER        (MAX_STA=5, NXSTEP=1, NZSTEP=1, NCHN=3) 
      !PARAMETER        (XSTART=-8000.0, XSTEP=500.0)
      !PARAMETER        (ZSTART=1000.0, ZSTEP=500.0)
      PARAMETER        (MAX_STA=5, NXSTEP=41, NZSTEP=21, NCHN=3) 
      PARAMETER        (XSTART=-10000.0, XSTEP=500.0)
      PARAMETER        (ZSTART=0.0, ZSTEP=500.0)
      PARAMETER        (VP=1500.0, VS=800.0, RHO=1700.0)

      DOUBLE PRECISION LAMBDA, MU, BETA, SXX, SYY, SZZ
      DOUBLE PRECISION PHI_M, THETA_M, RMS_M
      REAL             TIC, TOC
      INTEGER          NOBS, NSTA, ISTA, STAT, N1, IX, IY, IZ
      DOUBLE PRECISION RXX(MAX_STA), RYY(MAX_STA)
      DOUBLE PRECISION RZZ1(MAX_STA), RZZ2(MAX_STA)
      DOUBLE PRECISION BU(MAX_STA,2), LE(MAX_STA,2), LN(MAX_STA,2)
      DOUBLE PRECISION Y_O(MAX_STA*NCHN), DY_O(MAX_STA*NCHN)
      DOUBLE PRECISION X_M(4)
      DOUBLE PRECISION POOL(MAX_STA,20)
      CHARACTER*80     INFILE, STR_BETA

      INTEGER          IPID, NPID, IERR, SBLK, NX1, NX2

      !!!
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPID, IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, IPID, IERR)
      !!!
      SBLK = NXSTEP / NPID
      IF (SBLK * NPID < NXSTEP) THEN
        SBLK = SBLK + 1
      END IF
      NX1 = IPID * SBLK + 1
      NX2 = NX1 + SBLK - 1
      !PRINT *, IPID, NPID, NX1, NX2
      IF (NX2 > NXSTEP) THEN
        NX2 = NXSTEP
      END IF
      !!!

      MU = VS * VS * RHO
      LAMBDA = VP * VP * RHO  - 2 * MU 
      CALL GETARG(1, INFILE)
      CALL GETARG(2, STR_BETA)
      READ(STR_BETA, '(F15.8)') BETA
      !write(*,'(a,$)')' Please type the file name of input data: '
      !READ(*,*) INFILE
      !write(*,'(a,$)')' Please type the value of beta: '
      !READ(*,*) BETA
      OPEN(10,FILE=INFILE,STATUS='OLD')
      
      DO ISTA=1,MAX_STA
        READ(10,*,IOSTAT=STAT) RXX(ISTA), RYY(ISTA), 
     &   RZZ1(ISTA), RZZ2(ISTA),
     &   BU(ISTA,1), LE(ISTA,1), LN(ISTA,1),
     &   BU(ISTA,2), LE(ISTA,2), LN(ISTA,2)
        IF (STAT<0) THEN
          NSTA = ISTA-1
        END IF
      END DO
      CLOSE(10, STATUS='KEEP')
      
      NOBS = NSTA * NCHN

      DO ISTA=1,NSTA
        N1 = ISTA * 3
        Y_O(N1-2) = BU(ISTA,1)
        Y_O(N1-1) = LE(ISTA,1)
        Y_O(N1) = LN(ISTA,1)
        DY_O(N1-2) = BU(ISTA,2)
        DY_O(N1-1) = LE(ISTA,2)
        DY_O(N1) = LN(ISTA,2)
        !PRINT *,Y_O(N1-2),Y_O(N1-1),Y_O(N1) 
      END DO

      !DO IX=1,NXSTEP
      DO IX=NX1,NX2
        SXX = (IX-1) * XSTEP + XSTART
        DO IY=1,NXSTEP
          SYY = (IY-1) * XSTEP + XSTART
          DO IZ=1,NZSTEP
            SZZ = (IZ-1) * ZSTEP + ZSTART
            CALL CPU_TIME(TIC)
            IF (BETA>0) THEN
              CALL NLOKADA_3P(BETA, LAMBDA, MU, 
     *         NSTA, NCHN, NOBS,
     *         SXX, SYY, SZZ,
     *         RXX, RYY, RZZ1, RZZ2, Y_O, DY_O, 
     *         X_M, PHI_M, THETA_M, RMS_M)
            ELSE
              CALL NLOKADA_2P(LAMBDA, MU, 
     *         NSTA, NCHN, NOBS,
     *         SXX, SYY, SZZ, 
     *         RXX, RYY, RZZ1, RZZ2, Y_O, DY_O, 
     *         X_M, PHI_M, THETA_M, RMS_M)
            END IF
            !MISFIT(IX,IY,IZ) = RMS
            CALL CPU_TIME(TOC)
            WRITE(*,100) SXX, SYY, SZZ, RMS_M,
     *        PHI_M, THETA_M,
     *        X_M(1), X_M(2), X_M(3), X_M(4),
     *        TOC-TIC
            !CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
          END DO
        END DO
      END DO
      CALL MPI_FINALIZE(IERR)
      STOP
  100 format (3F12.1,F15.6,2F9.1,4E15.6,F12.4)
      END
