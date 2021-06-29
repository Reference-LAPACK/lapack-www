*
************************************************************************
*
      SUBROUTINE ECTPMV( UPLO, TRANS, DIAG, N, AP, X, INCX )
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      COMPLEX*16         X( * )
      COMPLEX            AP( * )
*     ..
*
*  Purpose
*  =======
*
*  ECTPMV performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x,   or   x := conjg( A' )*x,
*
*  where x is n element vector and A is an n by n unit, or non-unit,
*  upper or lower triangular matrix. Additional precision arithmetic is
*  used in the computation.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   x := A*x.
*
*              TRANS = 'T' or 't'   x := A'*x.
*
*              TRANS = 'C' or 'c'   x := conjg( A' )*x.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  AP     - COMPLEX          array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with  UPLO = 'U' or 'u', the array AP must
*           contain the upper triangular matrix packed sequentially,
*           column by column, so that AP( 1 ) contains a( 1, 1 ),
*           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
*           respectively, and so on.
*           Before entry with UPLO = 'L' or 'l', the array AP must
*           contain the lower triangular matrix packed sequentially,
*           column by column, so that AP( 1 ) contains a( 1, 1 ),
*           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
*           respectively, and so on.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x. At least double precision arithmetic is
*           used in the computation of x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 20-July-1986.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     .. Local Scalars ..
      INTEGER            I, INFO, IX, J, JX, K, KK, KX
      LOGICAL            NOCONJ, NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF      ( .NOT.LSAME( UPLO, 'U' ).AND.
     $          .NOT.LSAME( UPLO, 'L' )      ) THEN
         INFO = 1
      ELSE IF ( .NOT.LSAME( TRANS, 'N' ).AND.
     $          .NOT.LSAME( TRANS, 'T' ).AND.
     $          .NOT.LSAME( TRANS, 'C' )      ) THEN
         INFO = 2
      ELSE IF ( .NOT.LSAME( DIAG, 'U' ).AND.
     $          .NOT.LSAME( DIAG, 'N' )      ) THEN
         INFO = 3
      ELSE IF ( N.LT.0 ) THEN
         INFO = 4
      ELSE IF ( INCX.EQ.0 ) THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ECTPMV', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOCONJ = LSAME( TRANS, 'T' )
      NOUNIT = ( DIAG .EQ.'N' ).OR.( DIAG .EQ.'n' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of AP are
*     accessed sequentially with one pass through AP.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x:= A*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            K = 1
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + X( J )*AP( K )
                        K      = K      + 1
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*AP( K )
                     K = K + 1
                  ELSE
                     K = K + J
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     KK = K
                     DO 30, K = KK, KK + J - 2
                        X( IX ) = X( IX ) + X( JX )*AP( K )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*AP( K )
                     K = K + 1
                  ELSE
                     K = K + J
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            K = ( N*( N + 1 ) )/2
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + X( J )*AP( K )
                        K      = K      - 1
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*AP( K )
                     K = K - 1
                  ELSE
                     K = K - ( N - J + 1 )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     KK = K
                     DO 70, K = KK, KK - ( N - ( J + 1 ) ), -1
                        X( IX ) = X( IX ) + X( JX )*AP( K )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*AP( K )
                     K = K - 1
                  ELSE
                     K = K - ( N - J + 1 )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := A'*x  or  x := conjg( A' )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            K = ( N*( N + 1 ) )/2
            IF( INCX.EQ.1 )THEN
               DO 110, J = N, 1, -1
                  IF( NOCONJ )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )*AP( K )
                     K = K - 1
                     DO 90, I = J - 1, 1, -1
                        X( J ) = X( J ) + AP( K )*X( I )
                        K      = K      - 1
   90                CONTINUE
                  ELSE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*CONJG( AP( K ) )
                     K = K - 1
                     DO 100, I = J - 1, 1, -1
                        X( J ) = X( J ) + CONJG( AP( K ) )*X( I )
                        K      = K      - 1
  100                CONTINUE
                  END IF
  110          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 140, J = N, 1, -1
                  IX = JX
                  IF( NOCONJ )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*AP( K )
                     KK = K - 1
                     DO 120, K = KK, KK - J + 2, -1
                        IX      = IX      - INCX
                        X( JX ) = X( JX ) + AP( K )*X( IX )
  120                CONTINUE
                  ELSE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*CONJG( AP( K ) )
                     KK = K - 1
                     DO 130, K = KK, KK - J + 2, -1
                        IX      = IX      - INCX
                        X( JX ) = X( JX ) + CONJG( AP( K ) )*X( IX )
  130                CONTINUE
                  END IF
                  JX = JX - INCX
  140          CONTINUE
            END IF
         ELSE
            K = 1
            IF( INCX.EQ.1 )THEN
              DO 170, J = 1, N
                  IF( NOCONJ )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )*AP( K )
                     K = K + 1
                     DO 150, I = J + 1, N
                        X( J ) = X( J ) + AP( K )*X( I )
                        K      = K      + 1
  150                CONTINUE
                  ELSE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*CONJG( AP( K ) )
                     K = K + 1
                    DO 160, I = J + 1, N
                        X( J ) = X( J ) + CONJG( AP( K ) )*X( I )
                        K      = K      + 1
  160                CONTINUE
                  END IF
  170          CONTINUE
            ELSE
               JX = KX
               DO 200, J = 1, N
                  IX = JX
                  IF( NOCONJ )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*AP( K )
                     KK = K + 1
                     DO 180, K = KK, KK + N - ( J + 1 )
                        IX      = IX      + INCX
                        X( JX ) = X( JX ) + AP( K )*X( IX )
  180                CONTINUE
                  ELSE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*CONJG( AP( K ) )
                     KK = K + 1
                     DO 190, K = KK, KK + N - ( J + 1 )
                        IX      = IX      + INCX
                        X( JX ) = X( JX ) + CONJG( AP( K ) )*X( IX )
  190                CONTINUE
                  END IF
                  JX = JX + INCX
  200          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of ECTPMV.
*
      END
