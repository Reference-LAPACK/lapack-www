*
************************************************************************
*
      SUBROUTINE ECHPR ( UPLO, N, ALPHA, X, INCX, AP )
*     .. Scalar Arguments ..
      REAL               ALPHA
      INTEGER            INCX, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      COMPLEX*16         X( * )
      COMPLEX            AP( * )
*     ..
*
*  Purpose
*  =======
*
*  ECHPR   performs the hermitian rank 1 operation
*
*     A := alpha*x*conjg( x' ) + A,
*
*  where alpha is a real scalar, x is an n element vector and A is an
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
*  ALPHA  - REAL            .
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
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, J, JX, K, KK, KX
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, DBLE, REAL
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
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ECHPR ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.REAL( ZERO ) ) )
     $   RETURN
*
*     Set the start point in X if the increment is not unity.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
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
         IF( INCX.EQ.1 )THEN
            DO 20, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*DCONJG( X( J ) )
                  DO 10, I = 1, J - 1
                     AP( K ) = AP( K ) + X( I )*TEMP
                     K       = K       + 1
   10             CONTINUE
                  AP( K ) = REAL( AP( K ) ) + DBLE( X( J )*TEMP )
               ELSE
                  K       = K + J - 1
                  AP( K ) = REAL( AP( K ) )
               END IF
               K = K + 1
   20       CONTINUE
         ELSE
            JX = KX
            DO 40, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*DCONJG( X( JX ) )
                  IX   = KX
                  KK   = K
                  DO 30, K = KK, KK + J - 2
                     AP( K ) = AP( K ) + X( IX )*TEMP
                     IX      = IX      + INCX
   30             CONTINUE
                  AP( K ) = REAL( AP( K ) ) + DBLE( X( JX )*TEMP )
               ELSE
                  K       = K + J - 1
                  AP( K ) = REAL( AP( K ) )
               END IF
               K  = K  + 1
               JX = JX + INCX
   40       CONTINUE
         END IF
      ELSE
*
*        Form  A  when lower triangle is stored in AP.
*
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*DCONJG( X( J ) )
                  AP( K ) = REAL( AP( K ) ) + DBLE( TEMP*X( J ) )
                  K       = K               + 1
                  DO 50, I = J + 1, N
                     AP( K ) = AP( K ) + X( I )*TEMP
                     K       = K       + 1
   50             CONTINUE
               ELSE
                  AP( K ) = REAL( AP( K ) )
                  K       = K + N - J + 1
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP    = ALPHA*DCONJG( X( JX ) )
                  AP( K ) = REAL( AP( K ) ) + DBLE( TEMP*X( JX ) )
                  IX      = JX
                  KK      = K               + 1
                  DO 70, K = KK, KK + N - J - 1
                     IX      = IX      + INCX
                     AP( K ) = AP( K ) + X( IX )*TEMP
   70             CONTINUE
               ELSE
                  AP( K ) = REAL( AP( K ) )
                  K       = K + N - J + 1
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ECHPR .
*
      END
