      subroutine zblas2()
*
      CALL ZGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL ZGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL ZHEMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL ZHBMV ( UPLO, N, K, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
      CALL ZHPMV ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
*
      CALL ZTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*
      CALL ZTBMV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
*
      CALL ZTPMV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
*
      CALL ZTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*
      CALL ZTBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
*
      CALL ZTPSV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
*
      CALL ZGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*
      CALL ZGERC ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*
      CALL ZHER  ( UPLO, N, ALPHA, X, INCX, A, LDA )
*
      CALL ZHPR  ( UPLO, N, ALPHA, X, INCX, AP )
*
      CALL ZHER2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*
      CALL ZHPR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
*
      stop
      end
