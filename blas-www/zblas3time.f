      PROGRAM ZB3TIM
*
***********************************************************************
*
*  Timing program for the COMPLEX*16       Level 3 BLAS.
*
*  The program calls a REAL function SECOND with no arguments, which
*  is assumed to return the central-processor time in seconds from
*  some fixed starting-time.
*
*  The program must be driven by a short data file, which specifies
*  values for the matrix dimensions M, N and K, and for the leading
*  array dimension LDA (= LDB = LDC). The data file is read from unit
*  NIN. Results are written to unit NOUT. NIN and NOUT are set to 5 and
*  6 respectively in a PARAMETER statement, but can easily be changed.
*
*  The first 9 records of the data file are read using list-directed
*  input, the last 6 records are read using the format ( A ). An
*  annotated example of a data file can be obtained by deleting the
*  first 3 characters from the following 15 lines:
*  Data file for timing program for the COMPLEX*16       Level 3 BLAS
*  5                     Number of values of M (maximum 12)
*  32 64 128 256 512     The values of M
*  5                     Number of values of N (maximum 12)
*  32 64 128 256 512     The values of N
*  5                     Number of values of K (maximum 12)
*  32 64 128 256 512     The values of K
*  2                     Number of values of LDA=LDB=LDC (maximum 12)
*  512 513               The values of LDA=LDB=LDC
*  ZGEMM     T           ('T' or 't' means 'Time this routine')
*  ZSYMM     T
*  ZHEMM     T
*  ZSYRK     T
*  ZHERK     T
*  ZSYR2K    T
*  ZHER2K    T
*  ZTRMM     T
*  ZTRSM     T
*
*  The maximum value for M, N and K (NMAX) is defined in a PARAMETER
*  statement as 512. The maximum value of LDA is defined in a PARAMETER
*  statement as NMAX+32. Both can easily be changed.
*
*  Each routine is timed for all possible values of the options
*  TRANSA, TRANSB, TRANS, SIDE and UPLO. Option DIAG is always set
*  to 'N'. All elements of the matrices are initialised to non-zero
*  values. ALPHA and BETA are set to values other than 0, 1, -1.
*  LDB and LDC are always set to the same value as LDA. The routines
*  are timed for all combinations of applicable values of M, N, K
*  and LDA. Note that this can take a considerable amount of time.
*  It is important that the Level 3 BLAS should perform as well as
*  possible for all shapes and sizes of matrices.
*
***********************************************************************
*
*  -- Written on 1-July-1988.
*     Jeremy Du Croz and Mick Pont, NAG Central Office.
*
*     .. Parameters ..
      INTEGER            NIN, NOUT
      PARAMETER          ( NIN = 5, NOUT = 6 )
      INTEGER            MAXVAL, MXNLDA
      PARAMETER          ( MAXVAL = 12, MXNLDA = 12 )
      INTEGER            NMAX, LDAMAX
      PARAMETER          ( NMAX = 512, LDAMAX = NMAX + 32 )
      INTEGER            LA, LB, LC
      PARAMETER          ( LA = NMAX*LDAMAX, LB = LA, LC = LA )
      INTEGER            NTRANS, NSIDES, NUPLOS
      PARAMETER          ( NTRANS = 3, NSIDES = 2, NUPLOS = 2 )
      INTEGER            LDRES
      PARAMETER          ( LDRES = MAXVAL )
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 9 )
*     .. Local Scalars ..
      INTEGER            I, IKVAL, ILDA, ISIDE, ISUB, ITA, ITB, IUPLO,
     $                   LDA, MATCH, NKVALS, NLDA, NMVALS, NNVALS
      CHARACTER*80       LINE
*     .. Local Arrays ..
      COMPLEX*16         A( LA ), B( LB ), C( LC )
      DOUBLE PRECISION   RESLTS( LDRES, MAXVAL, MXNLDA )
      INTEGER            KVALS( MAXVAL ), LDAVAL( MXNLDA ),
     $                   MVALS( MAXVAL ), NVALS( MAXVAL )
      LOGICAL            TSTSUB( NSUBS )
      CHARACTER*1        SIDES( NSIDES ), TRANS( NTRANS ),
     $                   UPLOS( NUPLOS )
      CHARACTER*6        NAMES( NSUBS )
*     .. External Subroutines ..
      EXTERNAL           ZTIM01, ZTIM02, ZTIM03, ZTIM04, ZTIM05, ZTIM06,
     $                   ZTIM07, ZTIM08, ZTIM09, PRNTAB
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DBLE
*     .. Data statements ..
      DATA               NAMES/'ZGEMM ', 'ZSYMM ', 'ZHEMM ', 'ZSYRK ',
     $                   'ZHERK ', 'ZSYR2K', 'ZHER2K', 'ZTRMM ',
     $                   'ZTRSM '/
      DATA               TRANS/'N', 'T', 'C'/
      DATA               SIDES/'L', 'R'/
      DATA               UPLOS/'U', 'L'/
*     .. Executable Statements ..
*
*     Read driving file.
*
      READ( NIN, FMT = * )
      READ( NIN, FMT = * )NMVALS
      IF( NMVALS.GT.MAXVAL )THEN
         WRITE( NOUT, FMT = * )'Too many values of M: NMVALS = ', NMVALS
         STOP
      END IF
      READ( NIN, FMT = * )( MVALS( I ), I = 1, NMVALS )
      DO 10 I = 1, NMVALS
         IF( MVALS( I ).GT.NMAX )THEN
            WRITE( NOUT, FMT = * )'A value of M is too large: value is '
     $         , MVALS( I )
            STOP
         END IF
   10 CONTINUE
      READ( NIN, FMT = * )NNVALS
      IF( NNVALS.GT.MAXVAL )THEN
         WRITE( NOUT, FMT = * )'Too many values of N: NNVALS = ', NNVALS
         STOP
      END IF
      READ( NIN, FMT = * )( NVALS( I ), I = 1, NNVALS )
      DO 20 I = 1, NNVALS
         IF( NVALS( I ).GT.NMAX )THEN
            WRITE( NOUT, FMT = * )'A value of N is too large: value is '
     $         , NVALS( I )
            STOP
         END IF
   20 CONTINUE
      READ( NIN, FMT = * )NKVALS
      IF( NKVALS.GT.MAXVAL )THEN
         WRITE( NOUT, FMT = * )'Too many values of K: NKVALS = ', NKVALS
         STOP
      END IF
      READ( NIN, FMT = * )( KVALS( I ), I = 1, NKVALS )
      DO 30 I = 1, NKVALS
         IF( KVALS( I ).GT.NMAX )THEN
            WRITE( NOUT, FMT = * )'A value of K is too large: value is '
     $         , KVALS( I )
            STOP
         END IF
   30 CONTINUE
      READ( NIN, FMT = * )NLDA
      IF( NLDA.GT.MXNLDA )THEN
         WRITE( NOUT, FMT = * )'Too many values of LDA: NLDA = ', NLDA
         STOP
      END IF
      READ( NIN, FMT = * )( LDAVAL( I ), I = 1, NLDA )
      DO 40 I = 1, NLDA
         IF( LDAVAL( I ).GT.LDAMAX )THEN
            WRITE( NOUT, FMT = * )
     $         'A value of LDA is too large: value is ', LDAVAL( I )
            STOP
         END IF
   40 CONTINUE
*
*     Decide which routines are to be timed.
*
   50 READ( NIN, FMT = '(A)', END = 80 )LINE
      MATCH = 0
      DO 60 I = 1, NSUBS
         IF( NAMES( I ).EQ.LINE( 1: 6 ) )
     $      MATCH = I
   60 CONTINUE
      IF( MATCH.EQ.0 )THEN
         WRITE( NOUT, FMT = * )'Unknown subroutine name: ''',
     $      LINE( 1: 6 ), ''''
         STOP
      END IF
      I = 7
   70 IF( LINE( I: I ).EQ.' '.AND.I.LT.80 )THEN
         I = I + 1
         GO TO 70
      END IF
*
*     We will time routine if first non-blank character is 'T' or 't'.
*
      TSTSUB( MATCH ) = LINE( I: I ).EQ.'T'.OR.LINE( I: I ).EQ.'t'
      GO TO 50
   80 CONTINUE
*
*     Time each routine
*
      DO 430 ISUB = 1, NSUBS
         IF( TSTSUB( ISUB ) )THEN
*
*           Generate data
*
            DO 90 I = 1, LA
               A( I ) = DCMPLX( 1.0D0 + DBLE( I )/( LA + 1 ),
     $                  -1.0D0 + DBLE( I )/( LA + 1 ) )
   90       CONTINUE
            DO 100 I = 1, LB
               B( I ) = DCMPLX( 1.0D0 + DBLE( I )/( LB + 1 ),
     $                  -1.0D0 + DBLE( I )/( LB + 1 ) )
  100       CONTINUE
            DO 110 I = 1, LC
               C( I ) = DCMPLX( 1.0D0 + DBLE( I )/( LC + 1 ),
     $                  -1.0D0 + DBLE( I )/( LC + 1 ) )
  110       CONTINUE
*
*           Print header.
*
            WRITE( NOUT, FMT = 9999 )NAMES( ISUB )
            IF( NLDA.EQ.1 )THEN
               WRITE( NOUT, FMT = 9998 )LDAVAL( 1 )
            ELSE
               DO 120 I = 1, NLDA
                  WRITE( NOUT, FMT = 9997 )I, LDAVAL( I )
  120          CONTINUE
            END IF
            IF( NAMES( ISUB ).EQ.'ZGEMM ' )THEN
               DO 160 ITA = 1, NTRANS
                  DO 150 ITB = 1, NTRANS
                     DO 140 IKVAL = 1, NKVALS
                        DO 130 ILDA = 1, NLDA
                           LDA = LDAVAL( ILDA )
                           CALL ZTIM01( TRANS( ITA ), TRANS( ITB ),
     $                                  NMAX, MVALS, NMVALS, NVALS,
     $                                  NNVALS, KVALS( IKVAL ), A, LDA,
     $                                  B, LDA, C, LDA, LDRES,
     $                                  RESLTS( 1, 1, ILDA ) )
  130                   CONTINUE
                        IF( IKVAL.EQ.1 )
     $                     WRITE( NOUT, FMT = 9996 )TRANS( ITA ),
     $                     TRANS( ITB )
                        WRITE( NOUT, FMT = 9995 )KVALS( IKVAL )
                        CALL PRNTAB( 'M', 'N', MVALS, NMVALS, NVALS,
     $                               NNVALS, LDAVAL, NLDA, RESLTS,
     $                               LDRES, MAXVAL, NOUT )
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'ZSYMM ' )THEN
               DO 190 ISIDE = 1, NSIDES
                  DO 180 IUPLO = 1, NUPLOS
                     DO 170 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        CALL ZTIM02( SIDES( ISIDE ), UPLOS( IUPLO ),
     $                               NMAX, MVALS, NMVALS, NVALS, NNVALS,
     $                               A, LDA, B, LDA, C, LDA, LDRES,
     $                               RESLTS( 1, 1, ILDA ) )
  170                CONTINUE
                     WRITE( NOUT, FMT = 9994 )SIDES( ISIDE ),
     $                  UPLOS( IUPLO )
                     CALL PRNTAB( 'M', 'N', MVALS, NMVALS, NVALS,
     $                            NNVALS, LDAVAL, NLDA, RESLTS, LDRES,
     $                            MAXVAL, NOUT )
  180             CONTINUE
  190          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'ZHEMM ' )THEN
               DO 220 ISIDE = 1, NSIDES
                  DO 210 IUPLO = 1, NUPLOS
                     DO 200 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        CALL ZTIM03( SIDES( ISIDE ), UPLOS( IUPLO ),
     $                               NMAX, MVALS, NMVALS, NVALS, NNVALS,
     $                               A, LDA, B, LDA, C, LDA, LDRES,
     $                               RESLTS( 1, 1, ILDA ) )
  200                CONTINUE
                     WRITE( NOUT, FMT = 9993 )SIDES( ISIDE ),
     $                  UPLOS( IUPLO )
                     CALL PRNTAB( 'M', 'N', MVALS, NMVALS, NVALS,
     $                            NNVALS, LDAVAL, NLDA, RESLTS, LDRES,
     $                            MAXVAL, NOUT )
  210             CONTINUE
  220          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'ZSYRK ' )THEN
               DO 250 IUPLO = 1, NUPLOS
                  DO 240 ITB = 1, NTRANS
                     IF( TRANS( ITB ).NE.'C' )THEN
                        DO 230 ILDA = 1, NLDA
                           LDA = LDAVAL( ILDA )
                           CALL ZTIM04( UPLOS( IUPLO ), TRANS( ITB ),
     $                                  NMAX, NVALS, NNVALS, KVALS,
     $                                  NKVALS, A, LDA, C, LDA, LDRES,
     $                                  RESLTS( 1, 1, ILDA ) )
  230                   CONTINUE
                        WRITE( NOUT, FMT = 9992 )UPLOS( IUPLO ),
     $                     TRANS( ITB )
                        CALL PRNTAB( 'K', 'N', KVALS, NKVALS, NVALS,
     $                               NNVALS, LDAVAL, NLDA, RESLTS,
     $                               LDRES, MAXVAL, NOUT )
                     END IF
  240             CONTINUE
  250          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'ZHERK ' )THEN
               DO 280 IUPLO = 1, NUPLOS
                  DO 270 ITB = 1, NTRANS
                     IF( TRANS( ITB ).NE.'T' )THEN
                        DO 260 ILDA = 1, NLDA
                           LDA = LDAVAL( ILDA )
                           CALL ZTIM05( UPLOS( IUPLO ), TRANS( ITB ),
     $                                  NMAX, NVALS, NNVALS, KVALS,
     $                                  NKVALS, A, LDA, C, LDA, LDRES,
     $                                  RESLTS( 1, 1, ILDA ) )
  260                   CONTINUE
                        WRITE( NOUT, FMT = 9991 )UPLOS( IUPLO ),
     $                     TRANS( ITB )
                        CALL PRNTAB( 'K', 'N', KVALS, NKVALS, NVALS,
     $                               NNVALS, LDAVAL, NLDA, RESLTS,
     $                               LDRES, MAXVAL, NOUT )
                     END IF
  270             CONTINUE
  280          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'ZSYR2K' )THEN
               DO 310 IUPLO = 1, NUPLOS
                  DO 300 ITB = 1, NTRANS
                     IF( TRANS( ITB ).NE.'C' )THEN
                        DO 290 ILDA = 1, NLDA
                           LDA = LDAVAL( ILDA )
                           CALL ZTIM06( UPLOS( IUPLO ), TRANS( ITB ),
     $                                  NMAX, NVALS, NNVALS, KVALS,
     $                                  NKVALS, A, LDA, B, LDA, C, LDA,
     $                                  LDRES, RESLTS( 1, 1, ILDA ) )
  290                   CONTINUE
                        WRITE( NOUT, FMT = 9990 )UPLOS( IUPLO ),
     $                     TRANS( ITB )
                        CALL PRNTAB( 'K', 'N', KVALS, NKVALS, NVALS,
     $                               NNVALS, LDAVAL, NLDA, RESLTS,
     $                               LDRES, MAXVAL, NOUT )
                     END IF
  300             CONTINUE
  310          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'ZHER2K' )THEN
               DO 340 IUPLO = 1, NUPLOS
                  DO 330 ITB = 1, NTRANS
                     IF( TRANS( ITB ).NE.'T' )THEN
                        DO 320 ILDA = 1, NLDA
                           LDA = LDAVAL( ILDA )
                           CALL ZTIM07( UPLOS( IUPLO ), TRANS( ITB ),
     $                                  NMAX, NVALS, NNVALS, KVALS,
     $                                  NKVALS, A, LDA, B, LDA, C, LDA,
     $                                  LDRES, RESLTS( 1, 1, ILDA ) )
  320                   CONTINUE
                        WRITE( NOUT, FMT = 9989 )UPLOS( IUPLO ),
     $                     TRANS( ITB )
                        CALL PRNTAB( 'K', 'N', KVALS, NKVALS, NVALS,
     $                               NNVALS, LDAVAL, NLDA, RESLTS,
     $                               LDRES, MAXVAL, NOUT )
                     END IF
  330             CONTINUE
  340          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'ZTRMM ' )THEN
               DO 380 ISIDE = 1, NSIDES
                  DO 370 IUPLO = 1, NUPLOS
                     DO 360 ITA = 1, NTRANS
                        DO 350 ILDA = 1, NLDA
                           LDA = LDAVAL( ILDA )
                           CALL ZTIM08( SIDES( ISIDE ), UPLOS( IUPLO ),
     $                                  TRANS( ITA ), NMAX, MVALS,
     $                                  NMVALS, NVALS, NNVALS, A, LDA,
     $                                  B, LDA, LDRES,
     $                                  RESLTS( 1, 1, ILDA ) )
  350                   CONTINUE
                        WRITE( NOUT, FMT = 9988 )SIDES( ISIDE ),
     $                     UPLOS( IUPLO ), TRANS( ITA )
                        CALL PRNTAB( 'M', 'N', MVALS, NMVALS, NVALS,
     $                               NNVALS, LDAVAL, NLDA, RESLTS,
     $                               LDRES, MAXVAL, NOUT )
  360                CONTINUE
  370             CONTINUE
  380          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'ZTRSM' )THEN
               DO 420 ISIDE = 1, NSIDES
                  DO 410 IUPLO = 1, NUPLOS
                     DO 400 ITA = 1, NTRANS
                        DO 390 ILDA = 1, NLDA
                           LDA = LDAVAL( ILDA )
                           CALL ZTIM09( SIDES( ISIDE ), UPLOS( IUPLO ),
     $                                  TRANS( ITA ), NMAX, MVALS,
     $                                  NMVALS, NVALS, NNVALS, A, LDA,
     $                                  B, LDA, LDRES,
     $                                  RESLTS( 1, 1, ILDA ) )
  390                   CONTINUE
                        WRITE( NOUT, FMT = 9987 )SIDES( ISIDE ),
     $                     UPLOS( IUPLO ), TRANS( ITA )
                        CALL PRNTAB( 'M', 'N', MVALS, NMVALS, NVALS,
     $                               NNVALS, LDAVAL, NLDA, RESLTS,
     $                               LDRES, MAXVAL, NOUT )
  400                CONTINUE
  410             CONTINUE
  420          CONTINUE
            END IF
         END IF
         WRITE( NOUT, FMT = 9986 )
  430 CONTINUE
*
 9999 FORMAT( '1*** Speed of ', A, ' in megaflops ***' )
 9998 FORMAT( 5X, 'with LDA = ', I5 )
 9997 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
 9996 FORMAT( /1X, 'ZGEMM with TRANSA = ''', A, ''', TRANSB = ''', A,
     $      '''' )
 9995 FORMAT( /1X, 'K = ', I4, / )
 9994 FORMAT( /1X, 'ZSYMM with SIDE = ''', A, ''', UPLO = ''', A, '''',
     $      / )
 9993 FORMAT( /1X, 'ZHEMM with SIDE = ''', A, ''', UPLO = ''', A, '''',
     $      / )
 9992 FORMAT( /1X, 'ZSYRK with UPLO = ''', A, ''', TRANS = ''', A, '''',
     $      / )
 9991 FORMAT( /1X, 'ZHERK with UPLO = ''', A, ''', TRANS = ''', A, '''',
     $      / )
 9990 FORMAT( /1X, 'ZSYR2K with UPLO = ''', A, ''', TRANS = ''', A,
     $      '''', / )
 9989 FORMAT( /1X, 'ZHER2K with UPLO = ''', A, ''', TRANS = ''', A,
     $      '''', / )
 9988 FORMAT( /1X, 'ZTRMM with SIDE = ''', A, ''', UPLO = ''', A, ''',',
     $      ' TRANS = ''', A, '''', / )
 9987 FORMAT( /1X, 'ZTRSM with SIDE = ''', A, ''', UPLO = ''', A, ''',',
     $      ' TRANS = ''', A, '''', / )
 9986 FORMAT( ///// )
      END
*
      SUBROUTINE PRNTAB( LAB1, LAB2, M, LENM, N, LENN, LDAVAL, NLDA,
     $                   RESLTS, LDRES, SECRES, NOUT )
*     .. Scalar Arguments ..
      INTEGER            LDRES, LENM, LENN, NLDA, NOUT, SECRES
      CHARACTER*1        LAB1, LAB2
*     .. Array Arguments ..
      DOUBLE PRECISION   RESLTS( LDRES, SECRES, NLDA )
      INTEGER            LDAVAL( NLDA ), M( LENM ), N( LENN )
*     .. Local Scalars ..
      INTEGER            I, J, K
*     .. Executable Statements ..
      WRITE( NOUT, FMT = 9999 )LAB2, ( N( I ), I = 1, LENN )
      WRITE( NOUT, FMT = 9998 )LAB1
      DO 20 I = 1, LENM
         WRITE( NOUT, FMT = 9997 )M( I ),
     $      ( RESLTS( I, J, 1 ), J = 1, LENN )
         DO 10 K = 2, NLDA
            WRITE( NOUT, FMT = 9996 )( RESLTS( I, J, K ), J = 1, LENN )
   10    CONTINUE
         IF( NLDA.GT.1 )
     $      WRITE( NOUT, FMT = * )
   20 CONTINUE
      RETURN
*
 9999 FORMAT( 7X, A, I6, 11I9 )
 9998 FORMAT( 4X, A )
 9997 FORMAT( 1X, I4, 12F9.2 )
 9996 FORMAT( 5X, 12F9.2 )
      END
*
      SUBROUTINE ZTIM01( TRANSA, TRANSB, NMAX, MVALS, NMVALS, NVALS,
     $                   NNVALS, K, A, LDA, B, LDB, C, LDC, LDRES,
     $                   RESLTS )
*     Times ZGEMM.
*     .. Parameters ..
      COMPLEX*16         ALPHA, BETA
      PARAMETER          ( ALPHA = ( -0.9D0, 0.7D0 ),
     $                   BETA = ( 1.1D0, -0.3D0 ) )
*     .. Scalar Arguments ..
      INTEGER            K, LDA, LDB, LDC, LDRES, NMAX, NMVALS, NNVALS
      CHARACTER*1        TRANSA, TRANSB
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, NMAX ), B( LDB, NMAX ), C( LDC, NMAX )
      DOUBLE PRECISION   RESLTS( LDRES, NNVALS )
      INTEGER            MVALS( NMVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      DOUBLE PRECISION   S0, S1, S2, TIME
      INTEGER            I, IMVAL, INVAL, J, M, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           ZGEMM
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DBLE
*     .. Executable Statements ..
      DO 40 IMVAL = 1, NMVALS
         M = MVALS( IMVAL )
         DO 30 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise C each time to avoid overflow.
            DO 20 J = 1, N
               DO 10 I = 1, M
                  C( I, J ) = DCMPLX( 1.0D0 + DBLE( I + ( J - 1 )*N )/
     $                        ( M*N + 1 ), -1.0D0 +
     $                        DBLE( I + ( J - 1 )*N )/( M*N + 1 ) )
   10          CONTINUE
   20       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL ZGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                  BETA, C, LDC )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            NOP = 8*M*N*K
            IF( TIME.LE.0.0D0 )THEN
               RESLTS( IMVAL, INVAL ) = 0.0D0
            ELSE
               RESLTS( IMVAL, INVAL ) = NOP/( 1.0D6*TIME )
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
*
      SUBROUTINE ZTIM02( SIDE, UPLO, NMAX, MVALS, NMVALS, NVALS, NNVALS,
     $                   A, LDA, B, LDB, C, LDC, LDRES, RESLTS )
*     Times ZSYMM.
*     .. Parameters ..
      COMPLEX*16         ALPHA, BETA
      PARAMETER          ( ALPHA = ( -0.9D0, 0.7D0 ),
     $                   BETA = ( 1.1D0, -0.3D0 ) )
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDC, LDRES, NMAX, NMVALS, NNVALS
      CHARACTER*1        SIDE, UPLO
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, NMAX ), B( LDB, NMAX ), C( LDC, NMAX )
      DOUBLE PRECISION   RESLTS( LDRES, NNVALS )
      INTEGER            MVALS( NMVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      DOUBLE PRECISION   S0, S1, S2, TIME
      INTEGER            I, IMVAL, INVAL, J, M, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           ZSYMM
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DBLE
*     .. Executable Statements ..
      DO 40 IMVAL = 1, NMVALS
         M = MVALS( IMVAL )
         DO 30 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise C each time to avoid overflow.
            DO 20 J = 1, N
               DO 10 I = 1, M
                  C( I, J ) = DCMPLX( 1.0D0 + DBLE( I + ( J - 1 )*N )/
     $                        ( M*N + 1 ), -1.0D0 +
     $                        DBLE( I + ( J - 1 )*N )/( M*N + 1 ) )
   10          CONTINUE
   20       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL ZSYMM( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA,
     $                  C, LDC )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            IF( SIDE.EQ.'L' )THEN
               NOP = 8*M*M*N
            ELSE
               NOP = 8*M*N*N
            END IF
            IF( TIME.LE.0.0D0 )THEN
               RESLTS( IMVAL, INVAL ) = 0.0D0
            ELSE
               RESLTS( IMVAL, INVAL ) = NOP/( 1.0D6*TIME )
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
*
      SUBROUTINE ZTIM03( SIDE, UPLO, NMAX, MVALS, NMVALS, NVALS, NNVALS,
     $                   A, LDA, B, LDB, C, LDC, LDRES, RESLTS )
*     Times ZHEMM.
*     .. Parameters ..
      COMPLEX*16         ALPHA, BETA
      PARAMETER          ( ALPHA = ( -0.9D0, 0.7D0 ),
     $                   BETA = ( 1.1D0, -0.3D0 ) )
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDC, LDRES, NMAX, NMVALS, NNVALS
      CHARACTER*1        SIDE, UPLO
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, NMAX ), B( LDB, NMAX ), C( LDC, NMAX )
      DOUBLE PRECISION   RESLTS( LDRES, NNVALS )
      INTEGER            MVALS( NMVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      DOUBLE PRECISION   S0, S1, S2, TIME
      INTEGER            I, IMVAL, INVAL, J, M, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           ZHEMM
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DBLE
*     .. Executable Statements ..
      DO 40 IMVAL = 1, NMVALS
         M = MVALS( IMVAL )
         DO 30 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise C each time to avoid overflow.
            DO 20 J = 1, N
               DO 10 I = 1, M
                  C( I, J ) = DCMPLX( 1.0D0 + DBLE( I + ( J - 1 )*N )/
     $                        ( M*N + 1 ), -1.0D0 +
     $                        DBLE( I + ( J - 1 )*N )/( M*N + 1 ) )
   10          CONTINUE
   20       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL ZHEMM( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA,
     $                  C, LDC )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            IF( SIDE.EQ.'L' )THEN
               NOP = 8*M*M*N
            ELSE
               NOP = 8*M*N*N
            END IF
            IF( TIME.LE.0.0D0 )THEN
               RESLTS( IMVAL, INVAL ) = 0.0D0
            ELSE
               RESLTS( IMVAL, INVAL ) = NOP/( 1.0D6*TIME )
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
*
      SUBROUTINE ZTIM04( UPLO, TRANS, NMAX, NVALS, NNVALS, KVALS,
     $                   NKVALS, A, LDA, C, LDC, LDRES, RESLTS )
*     Times ZSYRK.
*     .. Parameters ..
      COMPLEX*16         ALPHA, BETA
      PARAMETER          ( ALPHA = ( -0.9D0, 0.7D0 ),
     $                   BETA = ( 1.1D0, -0.3D0 ) )
*     .. Scalar Arguments ..
      INTEGER            LDA, LDC, LDRES, NKVALS, NMAX, NNVALS
      CHARACTER*1        TRANS, UPLO
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, NMAX ), C( LDC, NMAX )
      DOUBLE PRECISION   RESLTS( LDRES, NNVALS )
      INTEGER            KVALS( NKVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      DOUBLE PRECISION   S0, S1, S2, TIME
      INTEGER            I, IKVAL, INVAL, J, K, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           ZSYRK
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DBLE
*     .. Executable Statements ..
      DO 40 IKVAL = 1, NKVALS
         K = KVALS( IKVAL )
         DO 30 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise C each time to avoid overflow.
            DO 20 J = 1, N
               DO 10 I = 1, N
                  C( I, J ) = DCMPLX( 1.0D0 + DBLE( I + ( J - 1 )*N )/
     $                        ( N*N + 1 ), -1.0D0 +
     $                        DBLE( I + ( J - 1 )*N )/( N*N + 1 ) )
   10          CONTINUE
   20       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL ZSYRK( UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            NOP = 4*N*N*K
            IF( TIME.LE.0.0D0 )THEN
               RESLTS( IKVAL, INVAL ) = 0.0D0
            ELSE
               RESLTS( IKVAL, INVAL ) = NOP/( 1.0D6*TIME )
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
*
      SUBROUTINE ZTIM05( UPLO, TRANS, NMAX, NVALS, NNVALS, KVALS,
     $                   NKVALS, A, LDA, C, LDC, LDRES, RESLTS )
*     Times ZHERK.
*     .. Parameters ..
      DOUBLE PRECISION   ALPHA, BETA
      PARAMETER          ( ALPHA = -0.9D0, BETA = 1.1D0 )
*     .. Scalar Arguments ..
      INTEGER            LDA, LDC, LDRES, NKVALS, NMAX, NNVALS
      CHARACTER*1        TRANS, UPLO
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, NMAX ), C( LDC, NMAX )
      DOUBLE PRECISION   RESLTS( LDRES, NNVALS )
      INTEGER            KVALS( NKVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      DOUBLE PRECISION   S0, S1, S2, TIME
      INTEGER            I, IKVAL, INVAL, J, K, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           ZHERK
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DBLE
*     .. Executable Statements ..
      DO 40 IKVAL = 1, NKVALS
         K = KVALS( IKVAL )
         DO 30 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise C each time to avoid overflow.
            DO 20 J = 1, N
               DO 10 I = 1, N
                  C( I, J ) = DCMPLX( 1.0D0 + DBLE( I + ( J - 1 )*N )/
     $                        ( N*N + 1 ), -1.0D0 +
     $                        DBLE( I + ( J - 1 )*N )/( N*N + 1 ) )
   10          CONTINUE
   20       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL ZHERK( UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            NOP = 4*N*N*K
            IF( TIME.LE.0.0D0 )THEN
               RESLTS( IKVAL, INVAL ) = 0.0D0
            ELSE
               RESLTS( IKVAL, INVAL ) = NOP/( 1.0D6*TIME )
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
*
      SUBROUTINE ZTIM06( UPLO, TRANS, NMAX, NVALS, NNVALS, KVALS,
     $                   NKVALS, A, LDA, B, LDB, C, LDC, LDRES, RESLTS )
*     Times ZSYR2K.
*     .. Parameters ..
      COMPLEX*16         ALPHA, BETA
      PARAMETER          ( ALPHA = ( -0.9D0, 0.7D0 ),
     $                   BETA = ( 1.1D0, -0.3D0 ) )
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDC, LDRES, NKVALS, NMAX, NNVALS
      CHARACTER*1        TRANS, UPLO
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, NMAX ), B( LDB, NMAX ), C( LDC, NMAX )
      DOUBLE PRECISION   RESLTS( LDRES, NNVALS )
      INTEGER            KVALS( NKVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      DOUBLE PRECISION   S0, S1, S2, TIME
      INTEGER            I, IKVAL, INVAL, J, K, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           ZSYR2K
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DBLE
*     .. Executable Statements ..
      DO 40 IKVAL = 1, NKVALS
         K = KVALS( IKVAL )
         DO 30 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise C each time to avoid overflow.
            DO 20 J = 1, N
               DO 10 I = 1, N
                  C( I, J ) = DCMPLX( 1.0D0 + DBLE( I + ( J - 1 )*N )/
     $                        ( N*N + 1 ), -1.0D0 +
     $                        DBLE( I + ( J - 1 )*N )/( N*N + 1 ) )
   10          CONTINUE
   20       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL ZSYR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA,
     $                   C, LDC )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            NOP = 8*N*N*K
            IF( TIME.LE.0.0D0 )THEN
               RESLTS( IKVAL, INVAL ) = 0.0D0
            ELSE
               RESLTS( IKVAL, INVAL ) = NOP/( 1.0D6*TIME )
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
*
      SUBROUTINE ZTIM07( UPLO, TRANS, NMAX, NVALS, NNVALS, KVALS,
     $                   NKVALS, A, LDA, B, LDB, C, LDC, LDRES, RESLTS )
*     Times ZHER2K.
*     .. Parameters ..
      COMPLEX*16         ALPHA
      PARAMETER          ( ALPHA = ( -0.9D0, 0.7D0 ) )
      DOUBLE PRECISION   BETA
      PARAMETER          ( BETA = 1.1D0 )
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDC, LDRES, NKVALS, NMAX, NNVALS
      CHARACTER*1        TRANS, UPLO
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, NMAX ), B( LDB, NMAX ), C( LDC, NMAX )
      DOUBLE PRECISION   RESLTS( LDRES, NNVALS )
      INTEGER            KVALS( NKVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      DOUBLE PRECISION   S0, S1, S2, TIME
      INTEGER            I, IKVAL, INVAL, J, K, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           ZHER2K
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DBLE
*     .. Executable Statements ..
      DO 40 IKVAL = 1, NKVALS
         K = KVALS( IKVAL )
         DO 30 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise C each time to avoid overflow.
            DO 20 J = 1, N
               DO 10 I = 1, N
                  C( I, J ) = DCMPLX( 1.0D0 + DBLE( I + ( J - 1 )*N )/
     $                        ( N*N + 1 ), -1.0D0 +
     $                        DBLE( I + ( J - 1 )*N )/( N*N + 1 ) )
   10          CONTINUE
   20       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL ZHER2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA,
     $                   C, LDC )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            NOP = 8*N*N*K
            IF( TIME.LE.0.0D0 )THEN
               RESLTS( IKVAL, INVAL ) = 0.0D0
            ELSE
               RESLTS( IKVAL, INVAL ) = NOP/( 1.0D6*TIME )
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
*
      SUBROUTINE ZTIM08( SIDE, UPLO, TRANSA, NMAX, MVALS, NMVALS, NVALS,
     $                   NNVALS, A, LDA, B, LDB, LDRES, RESLTS )
*     Times ZTRMM.
*     .. Parameters ..
      COMPLEX*16         ALPHA
      PARAMETER          ( ALPHA = ( -0.9D0, 0.7D0 ) )
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDRES, NMAX, NMVALS, NNVALS
      CHARACTER*1        SIDE, TRANSA, UPLO
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, NMAX ), B( LDB, NMAX )
      DOUBLE PRECISION   RESLTS( LDRES, NNVALS )
      INTEGER            MVALS( NMVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      DOUBLE PRECISION   S0, S1, S2, TIME
      INTEGER            I, IMVAL, INVAL, J, M, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           ZTRMM
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DBLE
*     .. Executable Statements ..
      DO 40 IMVAL = 1, NMVALS
         M = MVALS( IMVAL )
         DO 30 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise B each time to avoid overflow.
            DO 20 J = 1, N
               DO 10 I = 1, M
                  B( I, J ) = DCMPLX( 1.0D0 + DBLE( I + ( J - 1 )*N )/
     $                        ( M*N + 1 ), -1.0D0 +
     $                        DBLE( I + ( J - 1 )*N )/( M*N + 1 ) )
   10          CONTINUE
   20       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL ZTRMM( SIDE, UPLO, TRANSA, 'N', M, N, ALPHA, A, LDA, B,
     $                  LDB )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            IF( SIDE.EQ.'L' )THEN
               NOP = 4*M*M*N
            ELSE
               NOP = 4*M*N*N
            END IF
            IF( TIME.LE.0.0D0 )THEN
               RESLTS( IMVAL, INVAL ) = 0.0D0
            ELSE
               RESLTS( IMVAL, INVAL ) = NOP/( 1.0D6*TIME )
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
*
      SUBROUTINE ZTIM09( SIDE, UPLO, TRANSA, NMAX, MVALS, NMVALS, NVALS,
     $                   NNVALS, A, LDA, B, LDB, LDRES, RESLTS )
*     Times ZTRSM.
*     .. Parameters ..
      COMPLEX*16         ALPHA
      PARAMETER          ( ALPHA = ( -0.9D0, 0.7D0 ) )
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDRES, NMAX, NMVALS, NNVALS
      CHARACTER*1        SIDE, TRANSA, UPLO
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, NMAX ), B( LDB, NMAX )
      DOUBLE PRECISION   RESLTS( LDRES, NNVALS )
      INTEGER            MVALS( NMVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      DOUBLE PRECISION   S0, S1, S2, TIME
      INTEGER            I, IMVAL, INVAL, J, M, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           ZTRSM
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DBLE
*     .. Executable Statements ..
      DO 40 IMVAL = 1, NMVALS
         M = MVALS( IMVAL )
         DO 30 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise B each time to avoid overflow.
            DO 20 J = 1, N
               DO 10 I = 1, M
                  B( I, J ) = DCMPLX( 1.0D0 + DBLE( I + ( J - 1 )*N )/
     $                        ( M*N + 1 ), -1.0D0 +
     $                        DBLE( I + ( J - 1 )*N )/( M*N + 1 ) )
   10          CONTINUE
   20       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL ZTRSM( SIDE, UPLO, TRANSA, 'N', M, N, ALPHA, A, LDA, B,
     $                  LDB )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            IF( SIDE.EQ.'L' )THEN
               NOP = 4*M*M*N
            ELSE
               NOP = 4*M*N*N
            END IF
            IF( TIME.LE.0.0D0 )THEN
               RESLTS( IMVAL, INVAL ) = 0.0D0
            ELSE
               RESLTS( IMVAL, INVAL ) = NOP/( 1.0D6*TIME )
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END

