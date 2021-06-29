*
************************************************************************
*
      SUBROUTINE ECHPR2( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
*     .. Scalar Arguments ..
      COMPLEX            ALPHA
      INTEGER            INCX, INCY, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
      COMPLEX            AP( * )
*     ..
*
*  Purpose
*  =======
*
*  ECHPR2 performs the hermitian rank 2 operation
*
*     A := alpha*x*conjg( y' ) + conjg( alpha )*y*conjg( x' ) + A,
*
*  where alpha is a scalar, x and y are n element vectors and A is an
*  n by n hermitian matrix. Additional precision arithmetic is used in
*  the computation.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the matrix A is supplied in the packed
*           array AP as follows:
*
*              UPLO = 'U' or 'u'   The upper triangular part of A is
*                                  supplied in AP.
*
*              UPLO = 'L' or 'l'   The lower triangular part of A is
*                                  supplied in AP.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX         .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  AP     - COMPLEX          array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with  UPLO = 'U' or 'u', the array AP must
*           contain the upper triangular part of the hermitian matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
*           and a( 2, 2 ) respectively, and so on. On exit, the array
*           AP is overwritten by the upper triangular part of the
*           updated matrix.
*           Before entry with UPLO = 'L' or 'l', the array AP must
*           contain the lower triangular part of the hermitian matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
*           and a( 3, 1 ) respectively, and so on. On exit, the array
*           AP is overwritten by the lower triangular part of the
*           updated matrix.
*           Note that the imaginary parts of the diagonal elements need
*           not be set, they are assumed to be zero, and on exit they
*           are set to zero. At least double precision arithmetic is
*           used in the computation of A.
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
      COMPLEX*16         TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, DBLE, CMPLX, REAL
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF      ( .NOT.LSAME( UPLO, 'U' ).AND.
     $          .NOT.LSAME( UPLO, 'L' )      ) THEN
         INFO = 1
      ELSE IF ( N.LT.0 ) THEN
         INFO = 2
      ELSE IF ( INCX.EQ.0 ) THEN
         INFO = 5
      ELSE IF ( INCY.EQ.0 ) THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ECHPR2', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.CMPLX( ZERO ) ) )
     $   RETURN
*
*     Set up the start points in X and Y if the increments are not both
*     unity.
*
      IF( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
      END IF
*
*     Start the operations. In this version the elements of the array AP
*     are accessed sequentially with one pass through AP.
*
      K = 1
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  A  when upper triangle is stored in AP.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*DCONJG( Y( J ) )
                  TEMP2 = DCONJG( ALPHA*X( J ) )
                  DO 10, I = 1, J - 1
                     AP( K ) = AP( K ) + X( I )*TEMP1 + Y( I )*TEMP2
                     K       = K       + 1
   10             CONTINUE
                  AP( K ) = REAL( AP( K ) ) +
     $                      DBLE( X( J )*TEMP1 + Y( J )*TEMP2 )
               ELSE
                  K       = K + J - 1
                  AP( K ) = REAL( AP( K ) )
               END IF
               K = K + 1
   20       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 40, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*DCONJG( Y( JY ) )
                  TEMP2 = DCONJG( ALPHA*X( JX ) )
                  IX    = KX
                  IY    = KY
                  KK    = K
                  DO 30, K = KK, KK + J - 2
                     AP( K ) = AP( K ) + X( IX )*TEMP1 + Y( IY )*TEMP2
                     IX      = IX      + INCX
                     IY      = IY      + INCY
   30             CONTINUE
                  AP( K ) = REAL( AP( K ) ) +
     $                      DBLE( X( JX )*TEMP1 + Y( JY )*TEMP2 )
               ELSE
                  K       = K + J - 1
                  AP( K ) = REAL( AP( K ) )
               END IF
               K  = K  + 1
               JX = JX + INCX
               JY = JY + INCY
   40       CONTINUE
         END IF
      ELSE
*
*        Form  A  when lower triangle is stored in AP.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1   = ALPHA*DCONJG( Y( J ) )
                  TEMP2   = DCONJG( ALPHA*X( J ) )
                  AP( K ) = REAL( AP( K ) ) +
     $                      DBLE( X( J )*TEMP1 + Y( J )*TEMP2 )
                  K       = K               + 1
                  DO 50, I = J + 1, N
                     AP( K ) = AP( K ) + X( I )*TEMP1 + Y( I )*TEMP2
                     K       = K       + 1
   50             CONTINUE
               ELSE
                  AP( K ) = REAL( AP( K ) )
                  K       = K + N - J + 1
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1   = ALPHA*DCONJG( Y( JY ) )
                  TEMP2   = DCONJG( ALPHA*X( JX ) )
                  AP( K ) = REAL( AP( K ) ) +
     $                      DBLE( X( JX )*TEMP1 + Y( JY )*TEMP2 )
                  IX      = JX
                  IY      = JY
                  KK      = K               + 1
                  DO 70, K = KK, KK + N - J - 1
                     IX      = IX      + INCX
                     IY      = IY      + INCY
                     AP( K ) = AP( K ) + X( IX )*TEMP1 + Y( IY )*TEMP2
   70             CONTINUE
               ELSE
                  AP( K ) = REAL( AP( K ) )
                  K       = K + N - J + 1
               END IF
               JX = JX + INCX
               JY = JY + INCY
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ECHPR2.
*
      END
