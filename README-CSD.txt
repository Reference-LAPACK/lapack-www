== CS Decomposition (Brian Sutton, Randolph-Macon College, July 2010)

More info here: 

* link:http://faculty.rmc.edu/bsutton/csd-software.html[The fundamental CS decomposition]
* link:http://faculty.rmc.edu/bsutton/csd.html[Computing the complete CS decomposition]

=== SRC directory

.Contribution

- CUNCSD, DORCSD, SORCSD, ZUNCSD:
  Compute the CS decomposition of a block-partitioned orthogonal/unitary matrix.

- CUNBDB, DORBDB, SORBDB, ZUNBDB:
  Simultaneously bidiagonalizes the blocks of a partitioned orthogonal/unitary matrix.

- CBBCSD, DBBCSD, SBBCSD, ZBBCSD:
  Compute the CS decomposition of an orthogonal/unitary matrix in bidiagonal-block form.

- CLAPMR, DLAPMR, SLAMPR, ZLAPMR:
  Rearranges the rows of a matrix as specified by a permutation vector.

- DLARTGP, SLARTGP:
  Generate a plane rotation so that the "diagonal" (i.e., the scalar R) is nonnegative.

- DLARTGS, SLARTGS:
  Auxiliary subroutines.

.CUNCSD, DORCSD, SORCSD, ZUNCSD

LAPACK Driver Routines


The routine is recursive. It will re-enter only once if convenient.

[source,fortran]
----
      RECURSIVE SUBROUTINE DORCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS,
     $                             SIGNS, M, P, Q, X11, LDX11, X12,
     $                             LDX12, X21, LDX21, X22, LDX22, THETA,
     $                             U1, LDU1, U2, LDU2, V1T, LDV1T, V2T,
     $                             LDV2T, WORK, LWORK, IWORK, INFO )
*
*     Brian Sutton
*     Randolph-Macon College
*     July 2010
*
*     .. Scalar Arguments ..
      CHARACTER          JOBU1, JOBU2, JOBV1T, JOBV2T, SIGNS, TRANS
      INTEGER            INFO, LDU1, LDU2, LDV1T, LDV2T, LDX11, LDX12,
     $                   LDX21, LDX22, LRWORK, LWORK, M, P, Q
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   THETA( * )
      DOUBLE PRECISION   U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ),
     $                   V2T( LDV2T, * ), WORK( * ), X11( LDX11, * ),
     $                   X12( LDX12, * ), X21( LDX21, * ), X22( LDX22,
     $                   * )
*     ..
*
*  Purpose
*  =======
*
*  DORCSD computes the CS decomposition of an M-by-M partitioned
*  orthogonal matrix X:
*
*                                  [  I  0  0 |  0  0  0 ]
*                                  [  0  C  0 |  0 -S  0 ]
*      [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**T
*  X = [-----------] = [---------] [---------------------] [---------]   .
*      [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
*                                  [  0  S  0 |  0  C  0 ]
*                                  [  0  0  I |  0  0  0 ]
*
*  X11 is P-by-Q. The orthogonal matrices U1, U2, V1, and V2 are P-by-P,
*  (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are
*  R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in
*  which R = MIN(P,M-P,Q,M-Q).
*
*  Arguments
*  =========
*
*  JOBU1   (input) CHARACTER
*          = 'Y':      U1 is computed;
*          otherwise:  U1 is not computed.
*
*  JOBU2   (input) CHARACTER
*          = 'Y':      U2 is computed;
*          otherwise:  U2 is not computed.
*
*  JOBV1T  (input) CHARACTER
*          = 'Y':      V1T is computed;
*          otherwise:  V1T is not computed.
*
*  JOBV2T  (input) CHARACTER
*          = 'Y':      V2T is computed;
*          otherwise:  V2T is not computed.
*
*  TRANS   (input) CHARACTER
*          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major
*                      order;
*          otherwise:  X, U1, U2, V1T, and V2T are stored in column-
*                      major order.
*
*  SIGNS   (input) CHARACTER
*          = 'O':      The lower-left block is made nonpositive (the
*                      "other" convention);
*          otherwise:  The upper-right block is made nonpositive (the
*                      "default" convention).
*
*  M       (input) INTEGER
*          The number of rows and columns in X.
*
*  P       (input) INTEGER
*          The number of rows in X11 and X12. 0 <= P <= M.
*
*  Q       (input) INTEGER
*          The number of columns in X11 and X21. 0 <= Q <= M.
*
*  X       (input/workspace) DOUBLE PRECISION array, dimension (LDX,M)
*          On entry, the orthogonal matrix whose CSD is desired.
*
*  LDX     (input) INTEGER
*          The leading dimension of X. LDX >= MAX(1,M).
*
*  THETA   (output) DOUBLE PRECISION array, dimension (R), in which R =
*          MIN(P,M-P,Q,M-Q).
*          C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and
*          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ).
*
*  U1      (output) DOUBLE PRECISION array, dimension (P)
*          If JOBU1 = 'Y', U1 contains the P-by-P orthogonal matrix U1.
*
*  LDU1    (input) INTEGER
*          The leading dimension of U1. If JOBU1 = 'Y', LDU1 >=
*          MAX(1,P).
*
*  U2      (output) DOUBLE PRECISION array, dimension (M-P)
*          If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) orthogonal
*          matrix U2.
*
*  LDU2    (input) INTEGER
*          The leading dimension of U2. If JOBU2 = 'Y', LDU2 >=
*          MAX(1,M-P).
*
*  V1T     (output) DOUBLE PRECISION array, dimension (Q)
*          If JOBV1T = 'Y', V1T contains the Q-by-Q matrix orthogonal
*          matrix V1**T.
*
*  LDV1T   (input) INTEGER
*          The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >=
*          MAX(1,Q).
*
*  V2T     (output) DOUBLE PRECISION array, dimension (M-Q)
*          If JOBV2T = 'Y', V2T contains the (M-Q)-by-(M-Q) orthogonal
*          matrix V2**T.
*
*  LDV2T   (input) INTEGER
*          The leading dimension of V2T. If JOBV2T = 'Y', LDV2T >=
*          MAX(1,M-Q).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*          If INFO > 0 on exit, WORK(2:R) contains the values PHI(1),
*          ..., PHI(R-1) that, together with THETA(1), ..., THETA(R),
*          define the matrix in intermediate bidiagonal-block form
*          remaining after nonconvergence. INFO specifies the number
*          of nonzero PHI's.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the work array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace) INTEGER array, dimension (M-Q)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  DBBCSD did not converge. See the description of WORK
*                above for details.
*
*  Reference
*  =========
*
*  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.
*      Algorithms, 50(1):33-65, 2009.
*
*  ===================================================================
----
.CUNBDB, DORBDB, SORBDB, ZUNBDB

LAPACK Computational Routines

DORBDB simultaneously bidiagonalizes the blocks of an M-by-M partitioned
orthogonal matrix X.

[source,fortran]
----
      SUBROUTINE DORBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12,
     $                   X21, LDX21, X22, LDX22, THETA, PHI, TAUP1,
     $                   TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO )
*
*     Brian Sutton
*     Randolph-Macon College
*     July 2010
*
*     .. Scalar Arguments ..
      CHARACTER          SIGNS, TRANS
      INTEGER            INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P,
     $                   Q
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   PHI( * ), THETA( * )
      DOUBLE PRECISION   TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ),
     $                   WORK( * ), X11( LDX11, * ), X12( LDX12, * ),
     $                   X21( LDX21, * ), X22( LDX22, * )
*     ..
*
*  Purpose
*  =======
*
*  DORBDB simultaneously bidiagonalizes the blocks of an M-by-M
*  partitioned orthogonal matrix X:
*
*                                  [ B11 | B12 0  0 ]
*      [ X11 | X12 ]   [ P1 |    ] [  0  |  0 -I  0 ] [ Q1 |    ]**T
*  X = [-----------] = [---------] [----------------] [---------]   .
*      [ X21 | X22 ]   [    | P2 ] [ B21 | B22 0  0 ] [    | Q2 ]
*                                  [  0  |  0  0  I ]
*
*  X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is
*  not the case, then X must be transposed and/or permuted. This can be
*  done in constant time using the TRANS and SIGNS options. See DORCSD
*  for details.)
*
*  The orthogonal matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by-
*  (M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. They are
*  represented implicitly by Householder vectors.
*
*  B11, B12, B21, and B22 are Q-by-Q bidiagonal matrices represented
*  implicitly by angles THETA, PHI.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER
*          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major
*                      order;
*          otherwise:  X, U1, U2, V1T, and V2T are stored in column-
*                      major order.
*
*  SIGNS   (input) CHARACTER
*          = 'O':      The lower-left block is made nonpositive (the
*                      "other" convention);
*          otherwise:  The upper-right block is made nonpositive (the
*                      "default" convention).
*
*  M       (input) INTEGER
*          The number of rows and columns in X.
*
*  P       (input) INTEGER
*          The number of rows in X11 and X12. 0 <= P <= M.
*
*  Q       (input) INTEGER
*          The number of columns in X11 and X21. 0 <= Q <=
*          MIN(P,M-P,M-Q).
*
*  X11     (input/output) DOUBLE PRECISION array, dimension (LDX11,Q)
*          On entry, the top-left block of the orthogonal matrix to be
*          reduced. On exit, the form depends on TRANS:
*          If TRANS = 'N', then
*             the columns of tril(X11) specify reflectors for P1,
*             the rows of triu(X11,1) specify reflectors for Q1;
*          else TRANS = 'T', and
*             the rows of triu(X11) specify reflectors for P1,
*             the columns of tril(X11,-1) specify reflectors for Q1.
*
*  LDX11   (input) INTEGER
*          The leading dimension of X11. If TRANS = 'N', then LDX11 >=
*          P; else LDX11 >= Q.
*
*  X12     (input/output) DOUBLE PRECISION array, dimension (LDX12,M-Q)
*          On entry, the top-right block of the orthogonal matrix to
*          be reduced. On exit, the form depends on TRANS:
*          If TRANS = 'N', then
*             the rows of triu(X12) specify the first P reflectors for
*             Q2;
*          else TRANS = 'T', and
*             the columns of tril(X12) specify the first P reflectors
*             for Q2.
*
*  LDX12   (input) INTEGER
*          The leading dimension of X12. If TRANS = 'N', then LDX12 >=
*          P; else LDX11 >= M-Q.
*
*  X21     (input/output) DOUBLE PRECISION array, dimension (LDX21,Q)
*          On entry, the bottom-left block of the orthogonal matrix to
*          be reduced. On exit, the form depends on TRANS:
*          If TRANS = 'N', then
*             the columns of tril(X21) specify reflectors for P2;
*          else TRANS = 'T', and
*             the rows of triu(X21) specify reflectors for P2.
*
*  LDX21   (input) INTEGER
*          The leading dimension of X21. If TRANS = 'N', then LDX21 >=
*          M-P; else LDX21 >= Q.
*
*  X22     (input/output) DOUBLE PRECISION array, dimension (LDX22,M-Q)
*          On entry, the bottom-right block of the orthogonal matrix to
*          be reduced. On exit, the form depends on TRANS:
*          If TRANS = 'N', then
*             the rows of triu(X22(Q+1:M-P,P+1:M-Q)) specify the last
*             M-P-Q reflectors for Q2,
*          else TRANS = 'T', and
*             the columns of tril(X22(P+1:M-Q,Q+1:M-P)) specify the last
*             M-P-Q reflectors for P2.
*
*  LDX22   (input) INTEGER
*          The leading dimension of X22. If TRANS = 'N', then LDX22 >=
*          M-P; else LDX22 >= M-Q.
*
*  THETA   (output) DOUBLE PRECISION array, dimension (Q)
*          The entries of the bidiagonal blocks B11, B12, B21, B22 can
*          be computed from the angles THETA and PHI. See Further
*          Details.
*
*  PHI     (output) DOUBLE PRECISION array, dimension (Q-1)
*          The entries of the bidiagonal blocks B11, B12, B21, B22 can
*          be computed from the angles THETA and PHI. See Further
*          Details.
*
*  TAUP1   (output) DOUBLE PRECISION array, dimension (P)
*          The scalar factors of the elementary reflectors that define
*          P1.
*
*  TAUP2   (output) DOUBLE PRECISION array, dimension (M-P)
*          The scalar factors of the elementary reflectors that define
*          P2.
*
*  TAUQ1   (output) DOUBLE PRECISION array, dimension (Q)
*          The scalar factors of the elementary reflectors that define
*          Q1.
*
*  TAUQ2   (output) DOUBLE PRECISION array, dimension (M-Q)
*          The scalar factors of the elementary reflectors that define
*          Q2.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= M-Q.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The bidiagonal blocks B11, B12, B21, and B22 are represented
*  implicitly by angles THETA(1), ..., THETA(Q) and PHI(1), ...,
*  PHI(Q-1). B11 and B21 are upper bidiagonal, while B21 and B22 are
*  lower bidiagonal. Every entry in each bidiagonal band is a product
*  of a sine or cosine of a THETA with a sine or cosine of a PHI. See
*  [1] or DORCSD for details.
*
*  P1, P2, Q1, and Q2 are represented as products of elementary
*  reflectors. See DORCSD for details on generating P1, P2, Q1, and Q2
*  using DORGQR and DORGLQ.
*
*  Reference
*  =========
*
*  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.
*      Algorithms, 50(1):33-65, 2009.
*
*  ====================================================================
----

.CBBCSD, DBBCSD, SBBCSD, ZBBCSD

LAPACK Computational Routines

DBBCSD computes the CS decomposition of an orthogonal matrix in
bidiagonal-block form.

[source,fortran]
----
      SUBROUTINE DBBCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q,
     $                   THETA, PHI, U1, LDU1, U2, LDU2, V1T, LDV1T,
     $                   V2T, LDV2T, B11D, B11E, B12D, B12E, B21D, B21E,
     $                   B22D, B22E, WORK, LWORK, INFO )
*
*     Brian Sutton
*     Randolph-Macon College
*     July 2010
*
*     .. Scalar Arguments ..
      CHARACTER          JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS
      INTEGER            INFO, LDU1, LDU2, LDV1T, LDV2T, LWORK, M, P, Q
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B11D( * ), B11E( * ), B12D( * ), B12E( * ),
     $                   B21D( * ), B21E( * ), B22D( * ), B22E( * ),
     $                   PHI( * ), THETA( * ), WORK( * )
      DOUBLE PRECISION   U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ),
     $                   V2T( LDV2T, * )
*     ..
*
*  Purpose
*  =======
*
*  DBBCSD computes the CS decomposition of an orthogonal matrix in
*  bidiagonal-block form,
*
*
*      [ B11 | B12 0  0 ]
*      [  0  |  0 -I  0 ]
*  X = [----------------]
*      [ B21 | B22 0  0 ]
*      [  0  |  0  0  I ]
*
*                                [  C | -S  0  0 ]
*                    [ U1 |    ] [  0 |  0 -I  0 ] [ V1 |    ]**T
*                  = [---------] [---------------] [---------]   .
*                    [    | U2 ] [  S |  C  0  0 ] [    | V2 ]
*                                [  0 |  0  0  I ]
*
*  X is M-by-M, its top-left block is P-by-Q, and Q must be no larger
*  than P, M-P, or M-Q. (If Q is not the smallest index, then X must be
*  transposed and/or permuted. This can be done in constant time using
*  the TRANS and SIGNS options. See DORCSD for details.)
*
*  The bidiagonal matrices B11, B12, B21, and B22 are represented
*  implicitly by angles THETA(1:Q) and PHI(1:Q-1).
*
*  The orthogonal matrices U1, U2, V1T, and V2T are input/output.
*  The input matrices are pre- or post-multiplied by the appropriate
*  singular vector matrices.
*
*  Arguments
*  =========
*
*  JOBU1   (input) CHARACTER
*          = 'Y':      U1 is updated;
*          otherwise:  U1 is not updated.
*
*  JOBU2   (input) CHARACTER
*          = 'Y':      U2 is updated;
*          otherwise:  U2 is not updated.
*
*  JOBV1T  (input) CHARACTER
*          = 'Y':      V1T is updated;
*          otherwise:  V1T is not updated.
*
*  JOBV2T  (input) CHARACTER
*          = 'Y':      V2T is updated;
*          otherwise:  V2T is not updated.
*
*  TRANS   (input) CHARACTER
*          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major
*                      order;
*          otherwise:  X, U1, U2, V1T, and V2T are stored in column-
*                      major order.
*
*  M       (input) INTEGER
*          The number of rows and columns in X, the orthogonal matrix in
*          bidiagonal-block form.
*
*  P       (input) INTEGER
*          The number of rows in the top-left block of X. 0 <= P <= M.
*
*  Q       (input) INTEGER
*          The number of columns in the top-left block of X.
*          0 <= Q <= MIN(P,M-P,M-Q).
*
*  THETA   (input/output) DOUBLE PRECISION array, dimension (Q)
*          On entry, the angles THETA(1),...,THETA(Q) that, along with
*          PHI(1), ...,PHI(Q-1), define the matrix in bidiagonal-block
*          form. On exit, the angles whose cosines and sines define the
*          diagonal blocks in the CS decomposition.
*
*  PHI     (input/workspace) DOUBLE PRECISION array, dimension (Q-1)
*          The angles PHI(1),...,PHI(Q-1) that, along with THETA(1),...,
*          THETA(Q), define the matrix in bidiagonal-block form.
*
*  U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,P)
*          On entry, an LDU1-by-P matrix. On exit, U1 is postmultiplied
*          by the left singular vector matrix common to [ B11 ; 0 ] and
*          [ B12 0 0 ; 0 -I 0 0 ].
*
*  LDU1    (input) INTEGER
*          The leading dimension of the array U1.
*
*  U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,M-P)
*          On entry, an LDU2-by-(M-P) matrix. On exit, U2 is
*          postmultiplied by the left singular vector matrix common to
*          [ B21 ; 0 ] and [ B22 0 0 ; 0 0 I ].
*
*  LDU2    (input) INTEGER
*          The leading dimension of the array U2.
*
*  V1T     (input/output) DOUBLE PRECISION array, dimension (LDV1T,Q)
*          On entry, a LDV1T-by-Q matrix. On exit, V1T is premultiplied
*          by the transpose of the right singular vector
*          matrix common to [ B11 ; 0 ] and [ B21 ; 0 ].
*
*  LDV1T   (input) INTEGER
*          The leading dimension of the array V1T.
*
*  V2T     (input/output) DOUBLE PRECISION array, dimenison (LDV2T,M-Q)
*          On entry, a LDV2T-by-(M-Q) matrix. On exit, V2T is
*          premultiplied by the transpose of the right
*          singular vector matrix common to [ B12 0 0 ; 0 -I 0 ] and
*          [ B22 0 0 ; 0 0 I ].
*
*  LDV2T   (input) INTEGER
*          The leading dimension of the array V2T.
*
*  B11D    (output) DOUBLE PRECISION array, dimension (Q)
*          When DBBCSD converges, B11D contains the cosines of THETA(1),
*          ..., THETA(Q). If DBBCSD fails to converge, then B11D
*          contains the diagonal of the partially reduced top-left
*          block.
*
*  B11E    (output) DOUBLE PRECISION array, dimension (Q-1)
*          When DBBCSD converges, B11E contains zeros. If DBBCSD fails
*          to converge, then B11E contains the superdiagonal of the
*          partially reduced top-left block.
*
*  B12D    (output) DOUBLE PRECISION array, dimension (Q)
*          When DBBCSD converges, B12D contains the negative sines of
*          THETA(1), ..., THETA(Q). If DBBCSD fails to converge, then
*          B12D contains the diagonal of the partially reduced top-right
*          block.
*
*  B12E    (output) DOUBLE PRECISION array, dimension (Q-1)
*          When DBBCSD converges, B12E contains zeros. If DBBCSD fails
*          to converge, then B12E contains the subdiagonal of the
*          partially reduced top-right block.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= MAX(1,8*Q).
*
*          If LWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the WORK array,
*          returns this value as the first entry of the work array, and
*          no error message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if DBBCSD did not converge, INFO specifies the number
*                of nonzero entries in PHI, and B11D, B11E, etc.,
*                contain the partially reduced matrix.
*
*  Reference
*  =========
*
*  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.
*      Algorithms, 50(1):33-65, 2009.
*
*  Internal Parameters
*  ===================
*
*  TOLMUL  DOUBLE PRECISION, default = MAX(10,MIN(100,EPS**(-1/8)))
*          TOLMUL controls the convergence criterion of the QR loop.
*          Angles THETA(i), PHI(i) are rounded to 0 or PI/2 when they
*          are within TOLMUL*EPS of either bound.
*
*  ===================================================================
*
----

.CLAPMR, DLAPMR, SLAMPR, ZLAPMR

LAPACK Computational Routines

DLAPMR rearranges the rows of the M by N matrix X as specified by the
permutation K(1),K(2),...,K(M) of the integers 1,...,M. 

Routines similar to DLAPMT: existing DLAPMT works on columns, new DLAPMR on
rows.

[source,fortran]
----
      SUBROUTINE DLAPMR( FORWRD, M, N, X, LDX, K )
*
*     Originally DLAPMT
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     Adapted to DLAPMR by Brian Sutton
*     July 2010
*
*     .. Scalar Arguments ..
      LOGICAL            FORWRD
      INTEGER            LDX, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            K( * )
      DOUBLE PRECISION   X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  DLAPMR rearranges the rows of the M by N matrix X as specified
*  by the permutation K(1),K(2),...,K(M) of the integers 1,...,M.
*  If FORWRD = .TRUE.,  forward permutation:
*
*       X(K(I),*) is moved X(I,*) for I = 1,2,...,M.
*
*  If FORWRD = .FALSE., backward permutation:
*
*       X(I,*) is moved to X(K(I),*) for I = 1,2,...,M.
*
*  Arguments
*  =========
*
*  FORWRD  (input) LOGICAL
*          = .TRUE., forward permutation
*          = .FALSE., backward permutation
*
*  M       (input) INTEGER
*          The number of rows of the matrix X. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix X. N >= 0.
*
*  X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)
*          On entry, the M by N matrix X.
*          On exit, X contains the permuted matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X, LDX >= MAX(1,M).
*
*  K       (input/output) INTEGER array, dimension (M)
*          On entry, K contains the permutation vector. K is used as
*          internal workspace, but reset to its original value on
*          output.
*
*  =====================================================================
*
----

.DLARTGP, SLARTGP

LAPACK Computational Routines

DLARTGP generates a plane rotation (also known as Givens rotation). The
difference with existing DLARTG is that the sign is chosen so that the
"diagonal" (i.e., the scalar R) is nonnegative.  (Same difference as between
DLARFG and DLARFGP which use Householder.)

DLARTG/DLARTGP are slower, more accurate versions of the Level 1 BLAS routine
DROTG with some other slighter differences.


[source,fortran]
----
SUBROUTINE DLARTGP( F, G, CS, SN, R )
*
*     Originally DLARTG
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     Adapted to DLARTGP by Brian Sutton
*     July 2010
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
*     ..
*
*  Purpose
*  =======
*
*  DLARTGP generates a plane rotation so that
*
*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a slower, more accurate version of the BLAS1 routine DROTG,
*  with the following other differences:
*     F and G are unchanged on return.
*     If G=0, then CS=(+/-)1 and SN=0.
*     If F=0 and (G .ne. 0), then CS=0 and SN=(+/-)1.
*
*  The sign is chosen so that R >= 0.
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The first component of vector to be rotated.
*
*  G       (input) DOUBLE PRECISION
*          The second component of vector to be rotated.
*
*  CS      (output) DOUBLE PRECISION
*          The cosine of the rotation.
*
*  SN      (output) DOUBLE PRECISION
*          The sine of the rotation.
*
*  R       (output) DOUBLE PRECISION
*          The nonzero component of the rotated vector.
*
*  This version has a few statements commented out for thread safety
*  (machine parameters are computed on each entry). 10 feb 03, SJH.
*
*  =====================================================================
----

.CLARTGS, DLARTGS, SLARTGS, ZLARTGS

LAPACK Auxiliary Routines

[source,fortran]
----
     SUBROUTINE DLARTGS( X, Y, SIGMA, CS, SN )
      IMPLICIT NONE
*
*     Brian Sutton
*     Randolph-Macon College
*     July 2010
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION        CS, SIGMA, SN, X, Y
*     ..
*
*  Purpose
*  =======
*
*  DLARTGS generates a plane rotation designed to introduce a bulge in
*  Golub-Reinsch-style implicit QR iteration for the bidiagonal SVD
*  problem. X and Y are the top-row entries, and SIGMA is the shift.
*  The computed CS and SN define a plane rotation satisfying
*
*     [  CS  SN  ]  .  [ X^2 - SIGMA ]  =  [ R ],
*     [ -SN  CS  ]     [    X * Y    ]     [ 0 ]
*
*  with R nonnegative.  If X^2 - SIGMA and X * Y are 0, then the
*  rotation is by PI/2.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*          The (1,1) entry of an upper bidiagonal matrix.
*
*  Y       (input) DOUBLE PRECISION
*          The (1,2) entry of an upper bidiagonal matrix.
*
*  SIGMA   (input) DOUBLE PRECISION
*          The shift.
*
*  CS      (output) DOUBLE PRECISION
*          The cosine of the rotation.
*
*  SN      (output) DOUBLE PRECISION
*          The sine of the rotation.
*
*  ===================================================================
----
