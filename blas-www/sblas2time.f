      PROGRAM SB2TIM
*
***********************************************************************
*
*  Timing program for the REAL             Level 2 BLAS.
*
*  The program calls a REAL function SECOND with no arguments, which
*  is assumed to return the central-processor time in seconds from
*  some fixed starting-time.
*
*  The program must be driven by a short data file, which specifies
*  values for the matrix dimensions M and N, the bandwidth K, and
*  for the leading array dimension LDA. The data file is read from
*  unit NIN. Results are written to unit NOUT. NIN and NOUT are set
*  to 5 and 6 respectively in a PARAMETER statement, but can easily
*  be changed.
*
*  The first 9 records of the data file are read using list-directed
*  input, the last 16 records are read using the format ( A ). An
*  annotated example of a data file can be obtained by deleting the
*  first 3 characters from the following 27 lines:
*  Data file for timing program for the REAL             Level 2 BLAS
*  5                     Number of values of M (maximum 12)
*  32 64 128 256 512     The values of M
*  5                     Number of values of N (maximum 12)
*  32 64 128 256 512     The values of N
*  2                     Number of values of K (maximum 12)
*  8 32                  The values of K
*  2                     Number of values of INCX = INCY (maximum 4)
*  1 2                   The values of INCX = INCY
*  2                     Number of values of LDA (maximum 12)
*  512 513               The values of LDA
*  SGEMV     T           ('T' or 't' means 'Time this routine')
*  SGBMV     T
*  SSYMV     T
*  SSBMV     T
*  SSPMV     T
*  STRMV     T
*  STBMV     T
*  STPMV     T
*  STRSV     T
*  STBSV     T
*  STPSV     T
*  SGER      T
*  SSYR      T
*  SSPR      T
*  SSYR2     T
*  SSPR2     T
*
*  The maximum value for M and N (NMAX) is defined in a PARAMETER
*  statement as 512. The maximum value for INCX and INCY (INCMAX) is
*  defined in a PARAMETER statement as 32. The maximum value of LDA
*  is defined in a PARAMETER statement as NMAX+32. All three can
*  easily be changed. The maximum value for K is (NMAX-1)/2.
*
*  Each routine is timed for all possible values of the options
*  TRANS and UPLO. Option DIAG is always set to 'N'. All elements
*  of the matrices are initialised to non-zero values. ALPHA and
*  BETA are set to values other than 0, 1, -1. The routines are
*  timed for all combinations of applicable values of M, N, K, LDA
*  and INCX=INCY. Note that this can take a considerable amount of
*  time. It is important that the Level 2 BLAS should perform as well
*  as possible for all shapes and sizes of matrices.
*
***********************************************************************
*
*  -- Written on 28-August-1989.
*     Jeremy Du Croz and Mick Pont, NAG Central Office.
*
*     .. Parameters ..
      INTEGER            NIN, NOUT
      PARAMETER          ( NIN = 5, NOUT = 6 )
      INTEGER            MAXVAL, MXNINC, MXNLDA
      PARAMETER          ( MAXVAL = 12, MXNINC = 4, MXNLDA = 12 )
      INTEGER            NMAX, KMAX, INCMAX, LDAMAX
      PARAMETER          ( NMAX = 512, KMAX = ( NMAX - 1 )/2,
     $                   INCMAX = 32, LDAMAX = NMAX + 32 )
      INTEGER            LA, LX
      PARAMETER          ( LA = NMAX*LDAMAX, LX = NMAX*INCMAX )
      INTEGER            NTRANS, NUPLOS
      PARAMETER          ( NTRANS = 2, NUPLOS = 2 )
      INTEGER            LDRES
      PARAMETER          ( LDRES = MAXVAL )
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 16 )
*     .. Local Scalars ..
      INTEGER            I, IINC, ILDA, INCX, ISUB, ITA, IUPLO, J, LDA,
     $                   MATCH, NINCS, NKVALS, NLDA, NMVALS, NNVALS
      LOGICAL            IXANDY
      CHARACTER*80       LINE
*     .. Local Arrays ..
      REAL               A( LA ), RESLTS( LDRES, MAXVAL, MXNLDA,
     $                   MXNINC ), X( LX ), Y( LX )
      INTEGER            INCS( MXNINC ), KVALS( MAXVAL ),
     $                   LDAVAL( MXNLDA ), MVALS( MAXVAL ),
     $                   NVALS( MAXVAL )
      LOGICAL            TSTSUB( NSUBS )
      CHARACTER*1        TRANS( NTRANS ), UPLOS( NUPLOS )
      CHARACTER*6        NAMES( NSUBS )
*     .. External Subroutines ..
      EXTERNAL           PRNTAB, PRNTB1, STIM01, STIM02, STIM03, STIM04,
     $                   STIM05, STIM06, STIM07, STIM08, STIM09, STIM10,
     $                   STIM11, STIM12, STIM13, STIM14, STIM15, STIM16
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     .. Data statements ..
      DATA               NAMES/'SGEMV ', 'SGBMV ', 'SSYMV ', 'SSBMV ',
     $                   'SSPMV ', 'STRMV ', 'STBMV ', 'STPMV ',
     $                   'STRSV ', 'STBSV ', 'STPSV ', 'SGER  ',
     $                   'SSYR  ', 'SSPR  ', 'SSYR2 ', 'SSPR2 '/
      DATA               TRANS/'N', 'T'/
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
         IF( KVALS( I ).GT.KMAX )THEN
            WRITE( NOUT, FMT = * )'A value of K is too large: value is '
     $         , KVALS( I )
            STOP
         END IF
   30 CONTINUE
      READ( NIN, FMT = * )NINCS
      IF( NINCS.GT.MXNINC )THEN
         WRITE( NOUT, FMT = * )'Too many values of INCX: NINCS = ',
     $      NINCS
         STOP
      END IF
      READ( NIN, FMT = * )( INCS( I ), I = 1, NINCS )
      DO 40 I = 1, NINCS
         IF( INCS( I ).GT.INCMAX )THEN
            WRITE( NOUT, FMT = * )
     $         'A value of INCX is too large: value is ', INCS( I )
            STOP
         END IF
   40 CONTINUE
      READ( NIN, FMT = * )NLDA
      IF( NLDA.GT.MXNLDA )THEN
         WRITE( NOUT, FMT = * )'Too many values of LDA: NLDA = ', NLDA
         STOP
      END IF
      READ( NIN, FMT = * )( LDAVAL( I ), I = 1, NLDA )
      DO 50 I = 1, NLDA
         IF( LDAVAL( I ).GT.LDAMAX )THEN
            WRITE( NOUT, FMT = * )
     $         'A value of LDA is too large: value is ', LDAVAL( I )
            STOP
         END IF
   50 CONTINUE
*
*     Decide which routines are to be timed.
*
   60 READ( NIN, FMT = '(A)', END = 90 )LINE
      MATCH = 0
      DO 70 I = 1, NSUBS
         IF( NAMES( I ).EQ.LINE( 1: 6 ) )
     $      MATCH = I
   70 CONTINUE
      IF( MATCH.EQ.0 )THEN
         WRITE( NOUT, FMT = * )'Unknown subroutine name: ''',
     $      LINE( 1: 6 ), ''''
         STOP
      END IF
      I = 7
   80 IF( LINE( I: I ).EQ.' '.AND.I.LT.80 )THEN
         I = I + 1
         GO TO 80
      END IF
*
*     We will time routine if first non-blank character is 'T' or 't'.
*
      TSTSUB( MATCH ) = LINE( I: I ).EQ.'T'.OR.LINE( I: I ).EQ.'t'
      GO TO 60
   90 CONTINUE
*
*     Time each routine
*
      DO 630 ISUB = 1, NSUBS
         IF( TSTSUB( ISUB ) )THEN
*
*           Generate data
*
            DO 100 I = 1, LA
               A( I ) = 1.0E0 + REAL( I )/( LA + 1 )
  100       CONTINUE
            DO 110 I = 1, LX
               X( I ) = I*1.013E0
               Y( I ) = I*1.013E0
  110       CONTINUE
*
*           Print header.
*
            WRITE( NOUT, FMT = 9999 )NAMES( ISUB )
            IXANDY = ISUB.LE.5.OR.ISUB.EQ.12.OR.ISUB.EQ.15.OR.ISUB.EQ.16
            IF( NAMES( ISUB )( 3: 3 ).NE.'P' )THEN
               IF( NLDA*NINCS.EQ.1 )THEN
                  IF( IXANDY )THEN
                     WRITE( NOUT, FMT = 9998 )LDAVAL( 1 ), INCS( 1 )
                  ELSE
                     WRITE( NOUT, FMT = 9997 )LDAVAL( 1 ), INCS( 1 )
                  END IF
               ELSE
                  DO 130 I = 1, NLDA
                     DO 120 J = 1, NINCS
                        IF( IXANDY )THEN
                           WRITE( NOUT, FMT = 9994 )( I - 1 )*NINCS + J,
     $                        LDAVAL( I ), INCS( J )
                        ELSE
                           WRITE( NOUT, FMT = 9993 )( I - 1 )*NINCS + J,
     $                        LDAVAL( I ), INCS( J )
                        END IF
  120                CONTINUE
  130             CONTINUE
               END IF
            ELSE
               IF( NINCS.EQ.1 )THEN
                  IF( IXANDY )THEN
                     WRITE( NOUT, FMT = 9996 )INCS( 1 )
                  ELSE
                     WRITE( NOUT, FMT = 9995 )INCS( 1 )
                  END IF
               ELSE
                  DO 140 J = 1, NINCS
                     IF( IXANDY )THEN
                        WRITE( NOUT, FMT = 9992 )J, INCS( J )
                     ELSE
                        WRITE( NOUT, FMT = 9991 )J, INCS( J )
                     END IF
  140             CONTINUE
               END IF
            END IF
*
            IF( NAMES( ISUB ).EQ.'SGEMV ' )THEN
               DO 170 ITA = 1, NTRANS
                  DO 160 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 150 IINC = 1, NINCS
                        INCX = INCS( IINC )
                        CALL STIM01( TRANS( ITA ), NMAX, MVALS, NMVALS,
     $                               NVALS, NNVALS, A, LDA, INCX, X, Y,
     $                               LDRES, RESLTS( 1, 1, ILDA, IINC ) )
  150                CONTINUE
  160             CONTINUE
                  WRITE( NOUT, FMT = 9990 )TRANS( ITA )
                  CALL PRNTAB( 'M', 'N', MVALS, NMVALS, NVALS, NNVALS,
     $                         NINCS, NLDA, RESLTS, LDRES, MAXVAL,
     $                         MAXVAL, NOUT )
  170          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'SGBMV ' )THEN
               DO 200 ITA = 1, NTRANS
                  DO 190 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 180 IINC = 1, NINCS
                        INCX = INCS( IINC )
                        CALL STIM02( TRANS( ITA ), NMAX, KVALS, NKVALS,
     $                               NVALS, NNVALS, A, LDA, INCX, X, Y,
     $                               LDRES, RESLTS( 1, 1, ILDA, IINC ) )
  180                CONTINUE
  190             CONTINUE
                  WRITE( NOUT, FMT = 9989 )TRANS( ITA )
                  CALL PRNTAB( 'K', 'N', KVALS, NKVALS, NVALS, NNVALS,
     $                         NINCS, NLDA, RESLTS, LDRES, MAXVAL,
     $                         MAXVAL, NOUT )
  200          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'SSYMV ' )THEN
               DO 230 IUPLO = 1, NUPLOS
                  DO 220 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 210 IINC = 1, NINCS
                        INCX = INCS( IINC )
                        CALL STIM03( UPLOS( IUPLO ), NMAX, NVALS,
     $                               NNVALS, A, LDA, INCX, X, Y, LDRES,
     $                               RESLTS( 1, 1, ILDA, IINC ) )
  210                CONTINUE
  220             CONTINUE
                  WRITE( NOUT, FMT = 9988 )UPLOS( IUPLO )
                  CALL PRNTB1( 'N', NVALS, NNVALS, NINCS, NLDA, RESLTS,
     $                         LDRES, MAXVAL, MAXVAL, NOUT )
  230          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'SSBMV ' )THEN
               DO 260 IUPLO = 1, NUPLOS
                  DO 250 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 240 IINC = 1, NINCS
                        INCX = INCS( IINC )
                        CALL STIM04( UPLOS( IUPLO ), NMAX, KVALS,
     $                               NKVALS, NVALS, NNVALS, A, LDA,
     $                               INCX, X, Y, LDRES,
     $                               RESLTS( 1, 1, ILDA, IINC ) )
  240                CONTINUE
  250             CONTINUE
                  WRITE( NOUT, FMT = 9987 )UPLOS( IUPLO )
                  CALL PRNTAB( 'K', 'N', KVALS, NKVALS, NVALS, NNVALS,
     $                         NINCS, NLDA, RESLTS, LDRES, MAXVAL,
     $                         MAXVAL, NOUT )
  260          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'SSPMV ' )THEN
               DO 280 IUPLO = 1, NUPLOS
                  ILDA = 1
                  LDA = LDAVAL( ILDA )
                  DO 270 IINC = 1, NINCS
                     INCX = INCS( IINC )
                     CALL STIM05( UPLOS( IUPLO ), NMAX, NVALS, NNVALS,
     $                            A, LDA, INCX, X, Y, LDRES,
     $                            RESLTS( 1, 1, ILDA, IINC ) )
  270             CONTINUE
                  WRITE( NOUT, FMT = 9986 )UPLOS( IUPLO )
                  CALL PRNTB1( 'N', NVALS, NNVALS, NINCS, 1, RESLTS,
     $                         LDRES, MAXVAL, MAXVAL, NOUT )
  280          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'STRMV ' )THEN
               DO 320 IUPLO = 1, NUPLOS
                  DO 310 ITA = 1, NTRANS
                     DO 300 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        DO 290 IINC = 1, NINCS
                           INCX = INCS( IINC )
                           CALL STIM06( UPLOS( IUPLO ), TRANS( ITA ),
     $                                  NMAX, NVALS, NNVALS, A, LDA,
     $                                  INCX, X, LDRES,
     $                                  RESLTS( 1, 1, ILDA, IINC ) )
  290                   CONTINUE
  300                CONTINUE
                     WRITE( NOUT, FMT = 9985 )UPLOS( IUPLO ),
     $                  TRANS( ITA )
                     CALL PRNTB1( 'N', NVALS, NNVALS, NINCS, NLDA,
     $                            RESLTS, LDRES, MAXVAL, MAXVAL, NOUT )
  310             CONTINUE
  320          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'STBMV ' )THEN
               DO 360 IUPLO = 1, NUPLOS
                  DO 350 ITA = 1, NTRANS
                     DO 340 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        DO 330 IINC = 1, NINCS
                           INCX = INCS( IINC )
                           CALL STIM07( UPLOS( IUPLO ), TRANS( ITA ),
     $                                  NMAX, KVALS, NKVALS, NVALS,
     $                                  NNVALS, A, LDA, INCX, X, LDRES,
     $                                  RESLTS( 1, 1, ILDA, IINC ) )
  330                   CONTINUE
  340                CONTINUE
                     WRITE( NOUT, FMT = 9984 )UPLOS( IUPLO ),
     $                  TRANS( ITA )
                     CALL PRNTAB( 'K', 'N', KVALS, NKVALS, NVALS,
     $                            NNVALS, NINCS, NLDA, RESLTS, LDRES,
     $                            MAXVAL, MAXVAL, NOUT )
  350             CONTINUE
  360          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'STPMV ' )THEN
               DO 390 IUPLO = 1, NUPLOS
                  DO 380 ITA = 1, NTRANS
                     ILDA = 1
                     LDA = LDAVAL( ILDA )
                     DO 370 IINC = 1, NINCS
                        INCX = INCS( IINC )
                        CALL STIM08( UPLOS( IUPLO ), TRANS( ITA ), NMAX,
     $                               NVALS, NNVALS, A, LDA, INCX, X,
     $                               LDRES, RESLTS( 1, 1, ILDA, IINC ) )
  370                CONTINUE
                     WRITE( NOUT, FMT = 9983 )UPLOS( IUPLO ),
     $                  TRANS( ITA )
                     CALL PRNTB1( 'N', NVALS, NNVALS, NINCS, 1, RESLTS,
     $                            LDRES, MAXVAL, MAXVAL, NOUT )
  380             CONTINUE
  390          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'STRSV ' )THEN
               DO 430 IUPLO = 1, NUPLOS
                  DO 420 ITA = 1, NTRANS
                     DO 410 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        DO 400 IINC = 1, NINCS
                           INCX = INCS( IINC )
                           CALL STIM09( UPLOS( IUPLO ), TRANS( ITA ),
     $                                  NMAX, NVALS, NNVALS, A, LDA,
     $                                  INCX, X, LDRES,
     $                                  RESLTS( 1, 1, ILDA, IINC ) )
  400                   CONTINUE
  410                CONTINUE
                     WRITE( NOUT, FMT = 9982 )UPLOS( IUPLO ),
     $                  TRANS( ITA )
                     CALL PRNTB1( 'N', NVALS, NNVALS, NINCS, NLDA,
     $                            RESLTS, LDRES, MAXVAL, MAXVAL, NOUT )
  420             CONTINUE
  430          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'STBSV ' )THEN
               DO 470 IUPLO = 1, NUPLOS
                  DO 460 ITA = 1, NTRANS
                     DO 450 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        DO 440 IINC = 1, NINCS
                           INCX = INCS( IINC )
                           CALL STIM10( UPLOS( IUPLO ), TRANS( ITA ),
     $                                  NMAX, KVALS, NKVALS, NVALS,
     $                                  NNVALS, A, LDA, INCX, X, LDRES,
     $                                  RESLTS( 1, 1, ILDA, IINC ) )
  440                   CONTINUE
  450                CONTINUE
                     WRITE( NOUT, FMT = 9981 )UPLOS( IUPLO ),
     $                  TRANS( ITA )
                     CALL PRNTAB( 'K', 'N', KVALS, NKVALS, NVALS,
     $                            NNVALS, NINCS, NLDA, RESLTS, LDRES,
     $                            MAXVAL, MAXVAL, NOUT )
  460             CONTINUE
  470          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'STPSV ' )THEN
               DO 500 IUPLO = 1, NUPLOS
                  DO 490 ITA = 1, NTRANS
                     ILDA = 1
                     LDA = LDAVAL( ILDA )
                     DO 480 IINC = 1, NINCS
                        INCX = INCS( IINC )
                        CALL STIM11( UPLOS( IUPLO ), TRANS( ITA ), NMAX,
     $                               NVALS, NNVALS, A, LDA, INCX, X,
     $                               LDRES, RESLTS( 1, 1, ILDA, IINC ) )
  480                CONTINUE
                     WRITE( NOUT, FMT = 9980 )UPLOS( IUPLO ),
     $                  TRANS( ITA )
                     CALL PRNTB1( 'N', NVALS, NNVALS, NINCS, 1, RESLTS,
     $                            LDRES, MAXVAL, MAXVAL, NOUT )
  490             CONTINUE
  500          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'SGER  ' )THEN
               DO 520 ILDA = 1, NLDA
                  LDA = LDAVAL( ILDA )
                  DO 510 IINC = 1, NINCS
                     INCX = INCS( IINC )
                     CALL STIM12( NMAX, MVALS, NMVALS, NVALS, NNVALS, A,
     $                            LDA, INCX, X, Y, LDRES,
     $                            RESLTS( 1, 1, ILDA, IINC ) )
  510             CONTINUE
  520          CONTINUE
               WRITE( NOUT, FMT = 9979 )
               CALL PRNTAB( 'M', 'N', MVALS, NMVALS, NVALS, NNVALS,
     $                      NINCS, NLDA, RESLTS, LDRES, MAXVAL, MAXVAL,
     $                      NOUT )
            ELSE IF( NAMES( ISUB ).EQ.'SSYR  ' )THEN
               DO 550 IUPLO = 1, NUPLOS
                  DO 540 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 530 IINC = 1, NINCS
                        INCX = INCS( IINC )
                        CALL STIM13( UPLOS( IUPLO ), NMAX, NVALS,
     $                               NNVALS, A, LDA, INCX, X, LDRES,
     $                               RESLTS( 1, 1, ILDA, IINC ) )
  530                CONTINUE
  540             CONTINUE
                  WRITE( NOUT, FMT = 9978 )UPLOS( IUPLO )
                  CALL PRNTB1( 'N', NVALS, NNVALS, NINCS, NLDA, RESLTS,
     $                         LDRES, MAXVAL, MAXVAL, NOUT )
  550          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'SSPR  ' )THEN
               DO 570 IUPLO = 1, NUPLOS
                  ILDA = 1
                  LDA = LDAVAL( ILDA )
                  DO 560 IINC = 1, NINCS
                     INCX = INCS( IINC )
                     CALL STIM14( UPLOS( IUPLO ), NMAX, NVALS, NNVALS,
     $                            A, LDA, INCX, X, LDRES,
     $                            RESLTS( 1, 1, ILDA, IINC ) )
  560             CONTINUE
                  WRITE( NOUT, FMT = 9977 )UPLOS( IUPLO )
                  CALL PRNTB1( 'N', NVALS, NNVALS, NINCS, 1, RESLTS,
     $                         LDRES, MAXVAL, MAXVAL, NOUT )
  570          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'SSYR2 ' )THEN
               DO 600 IUPLO = 1, NUPLOS
                  DO 590 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 580 IINC = 1, NINCS
                        INCX = INCS( IINC )
                        CALL STIM15( UPLOS( IUPLO ), NMAX, NVALS,
     $                               NNVALS, A, LDA, INCX, X, Y, LDRES,
     $                               RESLTS( 1, 1, ILDA, IINC ) )
  580                CONTINUE
  590             CONTINUE
                  WRITE( NOUT, FMT = 9976 )UPLOS( IUPLO )
                  CALL PRNTB1( 'N', NVALS, NNVALS, NINCS, NLDA, RESLTS,
     $                         LDRES, MAXVAL, MAXVAL, NOUT )
  600          CONTINUE
            ELSE IF( NAMES( ISUB ).EQ.'SSPR2 ' )THEN
               DO 620 IUPLO = 1, NUPLOS
                  ILDA = 1
                  LDA = LDAVAL( ILDA )
                  DO 610 IINC = 1, NINCS
                     INCX = INCS( IINC )
                     CALL STIM16( UPLOS( IUPLO ), NMAX, NVALS, NNVALS,
     $                            A, LDA, INCX, X, Y, LDRES,
     $                            RESLTS( 1, 1, ILDA, IINC ) )
  610             CONTINUE
                  WRITE( NOUT, FMT = 9975 )UPLOS( IUPLO )
                  CALL PRNTB1( 'N', NVALS, NNVALS, NINCS, 1, RESLTS,
     $                         LDRES, MAXVAL, MAXVAL, NOUT )
  620          CONTINUE
            END IF
            WRITE( NOUT, FMT = 9974 )
         END IF
  630 CONTINUE
*
 9999 FORMAT( '1*** Speed of ', A, ' in megaflops ***' )
 9998 FORMAT( 5X, 'with LDA = ', I5, ' and INCX = INCY = ', I5 )
 9997 FORMAT( 5X, 'with LDA = ', I5, ' and INCX = ', I5 )
 9996 FORMAT( 5X, 'with INCX = INCY = ', I5 )
 9995 FORMAT( 5X, 'with INCX = ', I5 )
 9994 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5,
     $      ' and INCX = INCY = ', I5 )
 9993 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5, ' and INCX = ', I5 )
 9992 FORMAT( 5X, 'line ', I2, ' with INCX = INCY = ', I5 )
 9991 FORMAT( 5X, 'line ', I2, ' with INCX = ', I5 )
 9990 FORMAT( /1X, 'SGEMV with TRANS = ''', A, '''', / )
 9989 FORMAT( /1X, 'SGBMV with TRANS = ''', A, ''', M = N and KL = KU ',
     $      '= K', / )
 9988 FORMAT( /1X, 'SSYMV with UPLO = ''', A, '''', / )
 9987 FORMAT( /1X, 'SSBMV with UPLO = ''', A, '''', / )
 9986 FORMAT( /1X, 'SSPMV with UPLO = ''', A, '''', / )
 9985 FORMAT( /1X, 'STRMV with UPLO = ''', A, ''', TRANS = ''', A, '''',
     $      / )
 9984 FORMAT( /1X, 'STBMV with UPLO = ''', A, ''', TRANS = ''', A, '''',
     $      / )
 9983 FORMAT( /1X, 'STPMV with UPLO = ''', A, ''', TRANS = ''', A, '''',
     $      / )
 9982 FORMAT( /1X, 'STRSV with UPLO = ''', A, ''', TRANS = ''', A, '''',
     $      / )
 9981 FORMAT( /1X, 'STBSV with UPLO = ''', A, ''', TRANS = ''', A, '''',
     $      / )
 9980 FORMAT( /1X, 'STPSV with UPLO = ''', A, ''', TRANS = ''', A, '''',
     $      / )
 9979 FORMAT( /1X, 'SGER', / )
 9978 FORMAT( /1X, 'SSYR with UPLO = ''', A, '''', / )
 9977 FORMAT( /1X, 'SSPR with UPLO = ''', A, '''', / )
 9976 FORMAT( /1X, 'SSYR2 with UPLO = ''', A, '''', / )
 9975 FORMAT( /1X, 'SSPR2 with UPLO = ''', A, '''', / )
 9974 FORMAT( ///// )
      END
*
      SUBROUTINE PRNTAB( LAB1, LAB2, M, LENM, N, LENN, NINCS, NLDA,
     $                   RESLTS, LDRES, SECRES, THDRES, NOUT )
*     .. Scalar Arguments ..
      INTEGER            LDRES, LENM, LENN, NINCS, NLDA, NOUT, SECRES,
     $                   THDRES
      CHARACTER*1        LAB1, LAB2
*     .. Array Arguments ..
      REAL               RESLTS( LDRES, SECRES, THDRES, NINCS )
      INTEGER            M( LENM ), N( LENN )
*     .. Local Scalars ..
      INTEGER            I, J, K, L
*     .. Executable Statements ..
      WRITE( NOUT, FMT = 9999 )LAB2, ( N( I ), I = 1, LENN )
      WRITE( NOUT, FMT = 9998 )LAB1
      DO 30 I = 1, LENM
         WRITE( NOUT, FMT = 9997 )M( I ),
     $      ( RESLTS( I, J, 1, 1 ), J = 1, LENN )
         DO 20 K = 1, NLDA
            DO 10 L = 1, NINCS
               IF( K*L.GT.1 )
     $            WRITE( NOUT, FMT = 9996 )( RESLTS( I, J, K, L ),
     $            J = 1, LENN )
   10       CONTINUE
   20    CONTINUE
         IF( NLDA*NINCS.GT.1 )
     $      WRITE( NOUT, FMT = * )
   30 CONTINUE
      RETURN
*
 9999 FORMAT( 7X, A, I6, 11I9 )
 9998 FORMAT( 4X, A )
 9997 FORMAT( 1X, I4, 12F9.2 )
 9996 FORMAT( 5X, 12F9.2 )
      END
*
      SUBROUTINE PRNTB1( LAB1, N, LENN, NINCS, NLDA, RESLTS, LDRES,
     $                   SECRES, THDRES, NOUT )
*     .. Scalar Arguments ..
      INTEGER            LDRES, LENN, NINCS, NLDA, NOUT, SECRES, THDRES
      CHARACTER*1        LAB1
*     .. Array Arguments ..
      REAL               RESLTS( LDRES, SECRES, THDRES, NINCS )
      INTEGER            N( LENN )
*     .. Local Scalars ..
      INTEGER            I, J, K, L
*     .. Executable Statements ..
      WRITE( NOUT, FMT = 9999 )LAB1, ( N( I ), I = 1, LENN )
      WRITE( NOUT, FMT = * )
      DO 20 K = 1, NLDA
         DO 10 L = 1, NINCS
            WRITE( NOUT, FMT = 9998 )( RESLTS( 1, J, K, L ), J = 1,
     $         LENN )
   10    CONTINUE
   20 CONTINUE
      WRITE( NOUT, FMT = * )
      RETURN
*
 9999 FORMAT( 7X, A, I6, 11I9 )
 9998 FORMAT( 5X, 12F9.2 )
      END
*
      SUBROUTINE STIM01( TRANS, NMAX, MVALS, NMVALS, NVALS, NNVALS, A,
     $                   LDA, INCX, X, Y, LDRES, RESLTS )
*     Times SGEMV.
*     .. Parameters ..
      REAL               ALPHA, BETA
      PARAMETER          ( ALPHA = -0.9E0, BETA = 1.1E0 )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, LDRES, NMAX, NMVALS, NNVALS
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX ), Y( NMAX*INCX )
      INTEGER            MVALS( NMVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, IMVAL, INVAL, M, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           SGEMV
*     .. Executable Statements ..
      DO 30 IMVAL = 1, NMVALS
         M = MVALS( IMVAL )
         DO 20 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise Y each time to avoid overflow.
            DO 10 I = 1, N*INCX
               Y( I ) = I*1.013E0
   10       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL SGEMV( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y,
     $                  INCX )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            NOP = 2*M*N + 2*M
            IF( TIME.LE.0.0E0 )THEN
               RESLTS( IMVAL, INVAL ) = 0.0E0
            ELSE
               RESLTS( IMVAL, INVAL ) = NOP/( 1.0E6*TIME )
            END IF
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM02( TRANS, NMAX, KVALS, NKVALS, NVALS, NNVALS, A,
     $                   LDA, INCX, X, Y, LDRES, RESLTS )
*     Times SGBMV.
*     .. Parameters ..
      REAL               ALPHA, BETA
      PARAMETER          ( ALPHA = -0.9E0, BETA = 1.1E0 )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, LDRES, NKVALS, NMAX, NNVALS
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX ), Y( NMAX*INCX )
      INTEGER            KVALS( NKVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, IKVAL, INVAL, K, K1, M, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           SGBMV
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     .. Executable Statements ..
      DO 30 IKVAL = 1, NKVALS
         K = KVALS( IKVAL )
         DO 20 INVAL = 1, NNVALS
            N = NVALS( INVAL )
            M = N
*           Re-initialise Y each time to avoid overflow.
            DO 10 I = 1, N*INCX
               Y( I ) = I*1.013E0
   10       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL SGBMV( TRANS, M, N, K, K, ALPHA, A, LDA, X, INCX, BETA,
     $                  Y, INCX )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            K1 = MIN( K, N )
            NOP = 2*( ( 2*K1 + 1 )*N - K1*( K1 + 1 ) ) + 2*N
            IF( TIME.LE.0.0E0 )THEN
               RESLTS( IKVAL, INVAL ) = 0.0E0
            ELSE
               RESLTS( IKVAL, INVAL ) = NOP/( 1.0E6*TIME )
            END IF
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM03( UPLO, NMAX, NVALS, NNVALS, A, LDA, INCX, X, Y,
     $                   LDRES, RESLTS )
*     Times SSYMV.
*     .. Parameters ..
      REAL               ALPHA, BETA
      PARAMETER          ( ALPHA = -0.9E0, BETA = 1.1E0 )
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO
      INTEGER            INCX, LDA, LDRES, NMAX, NNVALS
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX ), Y( NMAX*INCX )
      INTEGER            NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, INVAL, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           SSYMV
*     .. Executable Statements ..
      DO 20 INVAL = 1, NNVALS
         N = NVALS( INVAL )
*           Re-initialise Y each time to avoid overflow.
         DO 10 I = 1, N*INCX
            Y( I ) = I*1.013E0
   10    CONTINUE
         S0 = SECOND( )
         S1 = SECOND( )
         CALL SSYMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCX )
         S2 = SECOND( )
         TIME = ( S2 - S1 ) - ( S1 - S0 )
         NOP = 2*N*N + 2*N
         IF( TIME.LE.0.0E0 )THEN
            RESLTS( 1, INVAL ) = 0.0E0
         ELSE
            RESLTS( 1, INVAL ) = NOP/( 1.0E6*TIME )
         END IF
   20 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM04( UPLO, NMAX, KVALS, NKVALS, NVALS, NNVALS, A,
     $                   LDA, INCX, X, Y, LDRES, RESLTS )
*     Times SSBMV.
*     .. Parameters ..
      REAL               ALPHA, BETA
      PARAMETER          ( ALPHA = -0.9E0, BETA = 1.1E0 )
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO
      INTEGER            INCX, LDA, LDRES, NKVALS, NMAX, NNVALS
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX ), Y( NMAX*INCX )
      INTEGER            KVALS( NKVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, IKVAL, INVAL, K, K1, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           SSBMV
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     .. Executable Statements ..
      DO 30 IKVAL = 1, NKVALS
         K = KVALS( IKVAL )
         DO 20 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise Y each time to avoid overflow.
            DO 10 I = 1, N*INCX
               Y( I ) = I*1.013E0
   10       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL SSBMV( UPLO, N, K, ALPHA, A, LDA, X, INCX, BETA, Y,
     $                  INCX )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            K1 = MIN( K, N )
            NOP = 2*( ( 2*K1 + 1 )*N - K1*( K1 + 1 ) ) + 2*N
            IF( TIME.LE.0.0E0 )THEN
               RESLTS( IKVAL, INVAL ) = 0.0E0
            ELSE
               RESLTS( IKVAL, INVAL ) = NOP/( 1.0E6*TIME )
            END IF
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM05( UPLO, NMAX, NVALS, NNVALS, A, LDA, INCX, X, Y,
     $                   LDRES, RESLTS )
*     Times SSPMV.
*     .. Parameters ..
      REAL               ALPHA, BETA
      PARAMETER          ( ALPHA = -0.9E0, BETA = 1.1E0 )
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO
      INTEGER            INCX, LDA, LDRES, NMAX, NNVALS
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX ), Y( NMAX*INCX )
      INTEGER            NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, INVAL, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           SSPMV
*     .. Executable Statements ..
      DO 20 INVAL = 1, NNVALS
         N = NVALS( INVAL )
*           Re-initialise Y each time to avoid overflow.
         DO 10 I = 1, N*INCX
            Y( I ) = I*1.013E0
   10    CONTINUE
         S0 = SECOND( )
         S1 = SECOND( )
         CALL SSPMV( UPLO, N, ALPHA, A, X, INCX, BETA, Y, INCX )
         S2 = SECOND( )
         TIME = ( S2 - S1 ) - ( S1 - S0 )
         NOP = 2*N*N + 2*N
         IF( TIME.LE.0.0E0 )THEN
            RESLTS( 1, INVAL ) = 0.0E0
         ELSE
            RESLTS( 1, INVAL ) = NOP/( 1.0E6*TIME )
         END IF
   20 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM06( UPLO, TRANS, NMAX, NVALS, NNVALS, A, LDA, INCX,
     $                   X, LDRES, RESLTS )
*     Times STRMV.
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, LDRES, NMAX, NNVALS
      CHARACTER*1        TRANS, UPLO
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX )
      INTEGER            NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, INVAL, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           STRMV
*     .. Executable Statements ..
      DO 20 INVAL = 1, NNVALS
         N = NVALS( INVAL )
*           Re-initialise X each time to avoid overflow.
         DO 10 I = 1, N*INCX
            X( I ) = I*1.013E0
   10    CONTINUE
         S0 = SECOND( )
         S1 = SECOND( )
         CALL STRMV( UPLO, TRANS, 'N', N, A, LDA, X, INCX )
         S2 = SECOND( )
         TIME = ( S2 - S1 ) - ( S1 - S0 )
         NOP = N*N
         IF( TIME.LE.0.0E0 )THEN
            RESLTS( 1, INVAL ) = 0.0E0
         ELSE
            RESLTS( 1, INVAL ) = NOP/( 1.0E6*TIME )
         END IF
   20 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM07( UPLO, TRANS, NMAX, KVALS, NKVALS, NVALS,
     $                   NNVALS, A, LDA, INCX, X, LDRES, RESLTS )
*     Times STBMV.
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, LDRES, NKVALS, NMAX, NNVALS
      CHARACTER*1        TRANS, UPLO
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX )
      INTEGER            KVALS( NKVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, IKVAL, INVAL, K, K1, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           STBMV
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     .. Executable Statements ..
      DO 30 IKVAL = 1, NKVALS
         K = KVALS( IKVAL )
         DO 20 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise X each time to avoid overflow.
            DO 10 I = 1, N*INCX
               X( I ) = I*1.013E0
   10       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL STBMV( UPLO, TRANS, 'N', N, K, A, LDA, X, INCX )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            K1 = MIN( K, N )
            NOP = ( 2*N*K1 - K1*( K1 + 1 ) ) + N
            IF( TIME.LE.0.0E0 )THEN
               RESLTS( IKVAL, INVAL ) = 0.0E0
            ELSE
               RESLTS( IKVAL, INVAL ) = NOP/( 1.0E6*TIME )
            END IF
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM08( UPLO, TRANS, NMAX, NVALS, NNVALS, A, LDA, INCX,
     $                   X, LDRES, RESLTS )
*     Times STPMV.
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, LDRES, NMAX, NNVALS
      CHARACTER*1        TRANS, UPLO
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX )
      INTEGER            NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, INVAL, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           STPMV
*     .. Executable Statements ..
      DO 20 INVAL = 1, NNVALS
         N = NVALS( INVAL )
*           Re-initialise X each time to avoid overflow.
         DO 10 I = 1, N*INCX
            X( I ) = I*1.013E0
   10    CONTINUE
         S0 = SECOND( )
         S1 = SECOND( )
         CALL STPMV( UPLO, TRANS, 'N', N, A, X, INCX )
         S2 = SECOND( )
         TIME = ( S2 - S1 ) - ( S1 - S0 )
         NOP = N*N
         IF( TIME.LE.0.0E0 )THEN
            RESLTS( 1, INVAL ) = 0.0E0
         ELSE
            RESLTS( 1, INVAL ) = NOP/( 1.0E6*TIME )
         END IF
   20 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM09( UPLO, TRANS, NMAX, NVALS, NNVALS, A, LDA, INCX,
     $                   X, LDRES, RESLTS )
*     Times STRSV.
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, LDRES, NMAX, NNVALS
      CHARACTER*1        TRANS, UPLO
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX )
      INTEGER            NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, INVAL, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           STRSV
*     .. Executable Statements ..
      DO 20 INVAL = 1, NNVALS
         N = NVALS( INVAL )
*           Re-initialise X each time to avoid overflow.
         DO 10 I = 1, N*INCX
            X( I ) = I*1.013E0
   10    CONTINUE
         S0 = SECOND( )
         S1 = SECOND( )
         CALL STRSV( UPLO, TRANS, 'N', N, A, LDA, X, INCX )
         S2 = SECOND( )
         TIME = ( S2 - S1 ) - ( S1 - S0 )
         NOP = N*N
         IF( TIME.LE.0.0E0 )THEN
            RESLTS( 1, INVAL ) = 0.0E0
         ELSE
            RESLTS( 1, INVAL ) = NOP/( 1.0E6*TIME )
         END IF
   20 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM10( UPLO, TRANS, NMAX, KVALS, NKVALS, NVALS,
     $                   NNVALS, A, LDA, INCX, X, LDRES, RESLTS )
*     Times STBSV.
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, LDRES, NKVALS, NMAX, NNVALS
      CHARACTER*1        TRANS, UPLO
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX )
      INTEGER            KVALS( NKVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, IKVAL, INVAL, K, K1, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           STBSV
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     .. Executable Statements ..
      DO 30 IKVAL = 1, NKVALS
         K = KVALS( IKVAL )
         DO 20 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise X each time to avoid overflow.
            DO 10 I = 1, N*INCX
               X( I ) = I*1.013E0
   10       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL STBSV( UPLO, TRANS, 'N', N, K, A, LDA, X, INCX )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            K1 = MIN( K, N )
            NOP = ( 2*N*K1 - K1*( K1 + 1 ) ) + N
            IF( TIME.LE.0.0E0 )THEN
               RESLTS( IKVAL, INVAL ) = 0.0E0
            ELSE
               RESLTS( IKVAL, INVAL ) = NOP/( 1.0E6*TIME )
            END IF
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM11( UPLO, TRANS, NMAX, NVALS, NNVALS, A, LDA, INCX,
     $                   X, LDRES, RESLTS )
*     Times STPSV.
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, LDRES, NMAX, NNVALS
      CHARACTER*1        TRANS, UPLO
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX )
      INTEGER            NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, INVAL, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           STPSV
*     .. Executable Statements ..
      DO 20 INVAL = 1, NNVALS
         N = NVALS( INVAL )
*           Re-initialise X each time to avoid overflow.
         DO 10 I = 1, N*INCX
            X( I ) = I*1.013E0
   10    CONTINUE
         S0 = SECOND( )
         S1 = SECOND( )
         CALL STPSV( UPLO, TRANS, 'N', N, A, X, INCX )
         S2 = SECOND( )
         TIME = ( S2 - S1 ) - ( S1 - S0 )
         NOP = N*N
         IF( TIME.LE.0.0E0 )THEN
            RESLTS( 1, INVAL ) = 0.0E0
         ELSE
            RESLTS( 1, INVAL ) = NOP/( 1.0E6*TIME )
         END IF
   20 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM12( NMAX, MVALS, NMVALS, NVALS, NNVALS, A, LDA,
     $                   INCX, X, Y, LDRES, RESLTS )
*     Times SGER.
*     .. Parameters ..
      REAL               ALPHA
      PARAMETER          ( ALPHA = -0.9E0 )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, LDRES, NMAX, NMVALS, NNVALS
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX ), Y( NMAX*INCX )
      INTEGER            MVALS( NMVALS ), NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, IMVAL, INVAL, J, M, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           SGER
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, REAL
*     .. Executable Statements ..
      DO 40 IMVAL = 1, NMVALS
         M = MVALS( IMVAL )
         DO 30 INVAL = 1, NNVALS
            N = NVALS( INVAL )
*           Re-initialise A each time to avoid overflow.
            DO 20 J = 1, N
               DO 10 I = 1, M
                  A( I, J ) = 1.0E0 + REAL( ( J - 1 )*N + I )/
     $                        REAL( M*N + 1 )
   10          CONTINUE
   20       CONTINUE
            S0 = SECOND( )
            S1 = SECOND( )
            CALL SGER( M, N, ALPHA, X, INCX, Y, INCX, A, LDA )
            S2 = SECOND( )
            TIME = ( S2 - S1 ) - ( S1 - S0 )
            NOP = 2*M*N + MIN( M, N )
            IF( TIME.LE.0.0E0 )THEN
               RESLTS( IMVAL, INVAL ) = 0.0E0
            ELSE
               RESLTS( IMVAL, INVAL ) = NOP/( 1.0E6*TIME )
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM13( UPLO, NMAX, NVALS, NNVALS, A, LDA, INCX, X,
     $                   LDRES, RESLTS )
*     Times SSYR.
*     .. Parameters ..
      REAL               ALPHA
      PARAMETER          ( ALPHA = -0.9E0 )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, LDRES, NMAX, NNVALS
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX )
      INTEGER            NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, INVAL, J, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           SSYR
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     .. Executable Statements ..
      DO 30 INVAL = 1, NNVALS
         N = NVALS( INVAL )
*           Re-initialise A each time to avoid overflow.
         DO 20 J = 1, N
            DO 10 I = 1, N
               A( I, J ) = 1.0E0 + REAL( ( J - 1 )*N + I )/
     $                     REAL( N*N + 1 )
   10       CONTINUE
   20    CONTINUE
         S0 = SECOND( )
         S1 = SECOND( )
         CALL SSYR( UPLO, N, ALPHA, X, INCX, A, LDA )
         S2 = SECOND( )
         TIME = ( S2 - S1 ) - ( S1 - S0 )
         NOP = N*( N + 1 ) + N
         IF( TIME.LE.0.0E0 )THEN
            RESLTS( 1, INVAL ) = 0.0E0
         ELSE
            RESLTS( 1, INVAL ) = NOP/( 1.0E6*TIME )
         END IF
   30 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM14( UPLO, NMAX, NVALS, NNVALS, A, LDA, INCX, X,
     $                   LDRES, RESLTS )
*     Times SSPR.
*     .. Parameters ..
      REAL               ALPHA
      PARAMETER          ( ALPHA = -0.9E0 )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, LDRES, NMAX, NNVALS
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX )
      INTEGER            NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, INVAL, J, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           SSPR
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     .. Executable Statements ..
      DO 30 INVAL = 1, NNVALS
         N = NVALS( INVAL )
*           Re-initialise A each time to avoid overflow.
         DO 20 J = 1, N
            DO 10 I = 1, N
               A( I, J ) = 1.0E0 + REAL( ( J - 1 )*N + I )/
     $                     REAL( N*N + 1 )
   10       CONTINUE
   20    CONTINUE
         S0 = SECOND( )
         S1 = SECOND( )
         CALL SSPR( UPLO, N, ALPHA, X, INCX, A )
         S2 = SECOND( )
         TIME = ( S2 - S1 ) - ( S1 - S0 )
         NOP = N*( N + 1 ) + N
         IF( TIME.LE.0.0E0 )THEN
            RESLTS( 1, INVAL ) = 0.0E0
         ELSE
            RESLTS( 1, INVAL ) = NOP/( 1.0E6*TIME )
         END IF
   30 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM15( UPLO, NMAX, NVALS, NNVALS, A, LDA, INCX, X, Y,
     $                   LDRES, RESLTS )
*     Times SSYR2.
*     .. Parameters ..
      REAL               ALPHA
      PARAMETER          ( ALPHA = -0.9E0 )
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO
      INTEGER            INCX, LDA, LDRES, NMAX, NNVALS
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX ), Y( NMAX*INCX )
      INTEGER            NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, INVAL, J, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           SSYR2
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     .. Executable Statements ..
      DO 30 INVAL = 1, NNVALS
         N = NVALS( INVAL )
*           Re-initialise A each time to avoid overflow.
         DO 20 J = 1, N
            DO 10 I = 1, N
               A( I, J ) = 1.0E0 + REAL( ( J - 1 )*N + I )/
     $                     REAL( N*N + 1 )
   10       CONTINUE
   20    CONTINUE
         S0 = SECOND( )
         S1 = SECOND( )
         CALL SSYR2( UPLO, N, ALPHA, X, INCX, Y, INCX, A, LDA )
         S2 = SECOND( )
         TIME = ( S2 - S1 ) - ( S1 - S0 )
         NOP = 2*( N + 1 )*N + 2*N
         IF( TIME.LE.0.0E0 )THEN
            RESLTS( 1, INVAL ) = 0.0E0
         ELSE
            RESLTS( 1, INVAL ) = NOP/( 1.0E6*TIME )
         END IF
   30 CONTINUE
      RETURN
      END
*
      SUBROUTINE STIM16( UPLO, NMAX, NVALS, NNVALS, A, LDA, INCX, X, Y,
     $                   LDRES, RESLTS )
*     Times SSPR2.
*     .. Parameters ..
      REAL               ALPHA
      PARAMETER          ( ALPHA = -0.9E0 )
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO
      INTEGER            INCX, LDA, LDRES, NMAX, NNVALS
*     .. Array Arguments ..
      REAL               A( LDA, NMAX ), RESLTS( LDRES, NNVALS ),
     $                   X( NMAX*INCX ), Y( NMAX*INCX )
      INTEGER            NVALS( NNVALS )
*     .. Local Scalars ..
      REAL               S0, S1, S2, TIME
      INTEGER            I, INVAL, J, N, NOP
*     .. External Functions ..
      REAL               SECOND
      EXTERNAL           SECOND
*     .. External Subroutines ..
      EXTERNAL           SSPR2
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     .. Executable Statements ..
      DO 30 INVAL = 1, NNVALS
         N = NVALS( INVAL )
*           Re-initialise A each time to avoid overflow.
         DO 20 J = 1, N
            DO 10 I = 1, N
               A( I, J ) = 1.0E0 + REAL( ( J - 1 )*N + I )/
     $                     REAL( N*N + 1 )
   10       CONTINUE
   20    CONTINUE
         S0 = SECOND( )
         S1 = SECOND( )
         CALL SSPR2( UPLO, N, ALPHA, X, INCX, Y, INCX, A )
         S2 = SECOND( )
         TIME = ( S2 - S1 ) - ( S1 - S0 )
         NOP = 2*( N + 1 )*N + 2*N
         IF( TIME.LE.0.0E0 )THEN
            RESLTS( 1, INVAL ) = 0.0E0
         ELSE
            RESLTS( 1, INVAL ) = NOP/( 1.0E6*TIME )
         END IF
   30 CONTINUE
      RETURN
      END
