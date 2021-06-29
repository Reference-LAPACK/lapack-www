      subroutine cblas2()
*
      CALL CGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL CGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL CHEMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL CHBMV ( UPLO, N, K, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL CHPMV ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
*
      CALL CTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*
      CALL CTBMV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
*
      CALL CTPMV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
*
      CALL CTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*
      CALL CTBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
*
      CALL CTPSV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
*
      CALL CGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*
      CALL CGERC ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*
      CALL CHER  ( UPLO, N, ALPHA, X, INCX, A, LDA )
*
      CALL CHPR  ( UPLO, N, ALPHA, X, INCX, AP )
*
      CALL CHER2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*
      CALL CHPR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
*
      stop
      end
