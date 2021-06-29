    MODULE blas_dense_fpinfo

      ! Fortran 95 Bindings Generic Interface Blocks
      ! .  for section 2.8.10 -- Environmental Enquiry
      ! .  using external procedures.
      ! Written on 14 December 2000
      ! .  Zohair Maany and Sven Hammarling, NAG Central Office.

      PRIVATE
      PUBLIC :: fpinfo

      INTERFACE fpinfo
        FUNCTION fpinfo_d(cmach,prec)
          USE blas_operator_arguments, ONLY : blas_cmach_type
          USE blas_precision, ONLY : wp => dp
          REAL (wp), INTENT (IN) :: prec
          TYPE (blas_cmach_type), INTENT (IN) :: cmach
          REAL (wp) :: fpinfo_d
        END FUNCTION fpinfo_d
        FUNCTION fpinfo_s(cmach,prec)
          USE blas_operator_arguments, ONLY : blas_cmach_type
          USE blas_precision, ONLY : wp => sp
          REAL (wp), INTENT (IN) :: prec
          TYPE (blas_cmach_type), INTENT (IN) :: cmach
          REAL (wp) :: fpinfo_s
        END FUNCTION fpinfo_s
      END INTERFACE
    END MODULE blas_dense_fpinfo
