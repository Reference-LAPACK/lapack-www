      SUBROUTINE DGELUB(M,N,A,LDA,IPIV,INFO,NB)
      DOUBLE PRECISION A(LDA,N)
      INTEGER IPIV(N)
*
*     This routine computes an LU factorization of an m-by-n
*     matrix A, using partial pivoting with row interchanges.
*
      INFO = 0
      DO 100, J = 1, N, Nb
         JB = MIN(N-J+1,NB)
*
*        Apply previous interchanges to current block.
*
         DO 20, I = 1, J-1
            IP = IPIV(I)
            IF (IP.NE.I) CALL DSWAP (JB,A(I,J),LDA,A(IP,J),LDA)
   20    CONTINUE
*
*        Compute superdiagonal block of U.
*
         CALL DTRSM ('Left','Lower','No transpose','Unit',
     $               J-1,JB,A(1,1),LDA,A(1,J),LDA)
*
*        Update diagonal and subdiagonal blocks
*
         CALL DGEMM ('No transpose','No transpose',
     $               M-J+1,JB,J-1,-1.0D0,A(J,1),LDA,
     $               A(1,J),LDA,1.0D0,A(J,J),LDA)
*
*        Factorize diagonal and subdiagonal blocks.
*
         CALL DGELU (M-J+1,JB,A(J,J),LDA,IPIV(J),INFO)
         DO 30, I = J, J+JB-1
            IPIV(I) = J - 1 + IPIV(I)
   30    CONTINUE
*
*        Test for singularity.
*
         IF (INFO.NE.0) GO TO 120
*
*        Apply interchanges to previous blocks.
*
         DO 40, I = J, J+JB-1
            IP = IPIV(J)
            IF (IP.NE.I) CALL DSWAP (J-1,A(I,1),LDA,A(IP,1),LDA)
   40    CONTINUE
  100 CONTINUE
      RETURN
*
  120 INFO = INFO + J - 1
      RETURN
      END
