      subroutine dblas2()
*
      CALL DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL DGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL DSBMV ( UPLO, N, K, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL DSPMV ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
*
      CALL DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*
      CALL DTBMV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
*
      CALL DTPMV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
*
      CALL DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*
      CALL DTBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
*
      CALL DTPSV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
*
      CALL DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*
      CALL DSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
*
      CALL DSPR  ( UPLO, N, ALPHA, X, INCX, AP )
*
      CALL DSYR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*
      CALL DSPR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
*
      stop
      end
