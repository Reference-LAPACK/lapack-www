      SUBROUTINE DGELU (M,N,A,LDA,IPIV,INFO)
      DOUBLE PRECISION A(LDA,N)
      INTEGER IPIV(N)
*
*     This routine computes an LU factorization of an m-by-n
*     matrix A, using partial pivoting with row interchanges.
*
      DOUBLE PRECISION T
*
      INFO = 0
      DO 100, J = 1, N
*
*        Apply previous interchanges to jth column.
*
         DO 20, I = 1, J-1
            IP = IPIV(I)
            T = A(I,J)
            A(I,J) = A(IP,J)
            A(IP,J) = T
   20    CONTINUE
*
*        Compute elements 1:j-1
*
         CALL DTRSV ('Lower','No transpose','Unit',J-1,
     $               A(1,1),LDA,A(1,J),1)
*
*        Update elements j:n
*
         CALL DGEMV ('No transpose',M-J+1,J-1,-1.0D0,A(J,1),LDA,
     $               A(1,J),1,1.0D0,A(J,J),1)
*
*        Find pivot.
*
         JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
         IPIV(J) = JP
         IF (A(JP,J).EQ.0.0D0) GO TO 120
*
*        Apply interchange to columns 1:j.
*
         IF (JP.NE.J) CALL DSWAP (J,A(J,1),LDA,A(JP,1),LDA)
*
*        Compute elements j+1:m.
*
         CALL DSCAL (M-J,1.0D0/A(J,J),A(J+1,J),1)
  100 CONTINUE
      RETURN
*
  120 INFO = J
      RETURN
      END
