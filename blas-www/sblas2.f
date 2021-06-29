      subroutine sblas2()
*
      CALL SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL SGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL SSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL SSBMV ( UPLO, N, K, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL SSPMV ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
*
      CALL STRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*
      CALL STBMV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
*
      CALL STPMV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
*
      CALL STRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*
      CALL STBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
*
      CALL STPSV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
*
      CALL SGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*
      CALL SSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
*
      CALL SSPR  ( UPLO, N, ALPHA, X, INCX, AP )
*
      CALL SSYR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*
      CALL SSPR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
*
      stop
      end
