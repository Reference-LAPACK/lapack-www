    MODULE blas_operator_arguments
      ! Module for defining BLAS operator arguments, cmach values,
      ! sparse matrix properties, and prec
      PRIVATE
      PUBLIC :: blas_trans_type, blas_no_trans, blas_trans, blas_conj_trans, &
       blas_uplo_type, blas_upper, blas_lower, &
       blas_diag_type, blas_non_unit_diag, blas_unit_diag, &
       blas_side_type, blas_left_side, blas_right_side, &
       blas_cmach_type, blas_base, blas_t, blas_rnd, blas_ieee, &
                blas_emin, blas_emax, blas_eps, blas_prec, blas_underflow, &
                blas_overflow, blas_sfmin, &
       blas_norm_type, blas_one_norm, blas_real_one_norm, &
                blas_two_norm, blas_frobenius_norm, blas_inf_norm, &
                blas_real_inf_norm, blas_max_norm, blas_real_max_norm, &
       blas_sort_type, blas_decreasing_order, blas_increasing_order, &
       blas_conj_type, blas_no_conj, blas_conj, &
       blas_jrot_type, blas_jrot_inner, blas_jrot_outer, &
                blas_jrot_sorted, &
       blas_prec_type, blas_single, blas_double, blas_indigenous, &
                blas_extra, &
       blas_symmetry_type, blas_general, blas_symmetric, &
                blas_hermitian, blas_triangular, blas_lower_triangular, &
                blas_upper_triangular, &
       blas_field_type, blas_complex, blas_real, blas_double_precision, &
                blas_single_precision, &
       blas_size_type, blas_num_rows, blas_num_cols, blas_num_nonzeros, &
       blas_valid_handle_type, blas_valid_handle, &
       blas_sparsity_optimization, blas_regular, blas_irregular, &
                blas_block, blas_unassembled
      ! .. Derived Type Declarations ..
      TYPE :: blas_trans_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_trans_type
      TYPE :: blas_uplo_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_uplo_type
      TYPE :: blas_diag_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_diag_type
      TYPE :: blas_side_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_side_type
      TYPE :: blas_cmach_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_cmach_type
      TYPE :: blas_norm_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_norm_type
      TYPE :: blas_sort_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_sort_type
      TYPE :: blas_conj_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_conj_type
      TYPE :: blas_jrot_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_jrot_type
      TYPE :: blas_prec_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_prec_type
      TYPE :: blas_symmetry_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_symmetry_type
      TYPE :: blas_field_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_field_type
      TYPE :: blas_size_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_size_type
      TYPE :: blas_valid_handle_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_valid_handle_type
      TYPE :: blas_sparsity_optimization_type
        PRIVATE
        INTEGER :: dummy
      END TYPE blas_sparsity_optimization_type
      ! .. Dependents ..
      TYPE (blas_trans_type), PARAMETER :: blas_no_trans = blas_trans_type(1)
      TYPE (blas_trans_type), PARAMETER :: blas_trans = blas_trans_type(2)
      TYPE (blas_trans_type), PARAMETER :: blas_conj_trans = blas_trans_type(3)
      TYPE (blas_cmach_type), PARAMETER :: blas_base = blas_cmach_type(1)
      TYPE (blas_uplo_type), PARAMETER :: blas_upper = blas_uplo_type(1)
      TYPE (blas_uplo_type), PARAMETER :: blas_lower= blas_uplo_type(2)
      TYPE (blas_diag_type), PARAMETER :: blas_non_unit_diag = blas_diag_type(1)
      TYPE (blas_diag_type), PARAMETER :: blas_unit_diag = blas_diag_type(2)
      TYPE (blas_side_type), PARAMETER :: blas_left_side = blas_side_type(1)
      TYPE (blas_side_type), PARAMETER :: blas_right_side = blas_side_type(2)
      TYPE (blas_cmach_type), PARAMETER :: blas_t = blas_cmach_type(2)
      TYPE (blas_cmach_type), PARAMETER :: blas_rnd = blas_cmach_type(3)
      TYPE (blas_cmach_type), PARAMETER :: blas_ieee = blas_cmach_type(4)
      TYPE (blas_cmach_type), PARAMETER :: blas_emin = blas_cmach_type(5)
      TYPE (blas_cmach_type), PARAMETER :: blas_emax = blas_cmach_type(6)
      TYPE (blas_cmach_type), PARAMETER :: blas_eps = blas_cmach_type(7)
      TYPE (blas_cmach_type), PARAMETER :: blas_prec = blas_cmach_type(8)
      TYPE (blas_cmach_type), PARAMETER :: blas_underflow = blas_cmach_type(9)
      TYPE (blas_cmach_type), PARAMETER :: blas_overflow = blas_cmach_type(10)
      TYPE (blas_cmach_type), PARAMETER :: blas_sfmin = blas_cmach_type(11)
      TYPE (blas_norm_type), PARAMETER :: blas_one_norm = blas_norm_type(1)
      TYPE (blas_norm_type), PARAMETER :: blas_real_one_norm = blas_norm_type(2)
      TYPE (blas_norm_type), PARAMETER :: blas_two_norm = blas_norm_type(3)
      TYPE (blas_norm_type), PARAMETER ::  &
           blas_frobenius_norm = blas_norm_type(4)
      TYPE (blas_norm_type), PARAMETER :: blas_inf_norm = blas_norm_type(5)
      TYPE (blas_norm_type), PARAMETER :: blas_real_inf_norm = blas_norm_type(6)
      TYPE (blas_norm_type), PARAMETER :: blas_max_norm = blas_norm_type(7)
      TYPE (blas_norm_type), PARAMETER :: blas_real_max_norm = blas_norm_type(8)
      TYPE (blas_sort_type), PARAMETER ::  &
           blas_decreasing_order = blas_sort_type(1)
      TYPE (blas_sort_type), PARAMETER ::  &
           blas_increasing_order = blas_sort_type(2)
      TYPE (blas_conj_type), PARAMETER :: blas_conj = blas_conj_type(1)
      TYPE (blas_conj_type), PARAMETER :: blas_no_conj = blas_conj_type(2)
      TYPE (blas_jrot_type), PARAMETER :: blas_jrot_inner = blas_jrot_type(1)
      TYPE (blas_jrot_type), PARAMETER :: blas_jrot_outer = blas_jrot_type(2)
      TYPE (blas_jrot_type), PARAMETER :: blas_jrot_sorted = blas_jrot_type(3)
      TYPE (blas_prec_type), PARAMETER :: blas_single = blas_prec_type(1)
      TYPE (blas_prec_type), PARAMETER :: blas_double = blas_prec_type(2)
      TYPE (blas_prec_type), PARAMETER :: blas_indigenous = blas_prec_type(3)
      TYPE (blas_prec_type), PARAMETER :: blas_extra = blas_prec_type(4)
      TYPE (blas_symmetry_type), PARAMETER ::  &
           blas_general = blas_symmetry_type(1)
      TYPE (blas_symmetry_type), PARAMETER ::  &
           blas_symmetric = blas_symmetry_type(2)
      TYPE (blas_symmetry_type), PARAMETER ::  &
           blas_hermitian = blas_symmetry_type(3)
      TYPE (blas_symmetry_type), PARAMETER ::  &
           blas_triangular = blas_symmetry_type(4)
      TYPE (blas_symmetry_type), PARAMETER ::  &
           blas_lower_triangular = blas_symmetry_type(5)
      TYPE (blas_symmetry_type), PARAMETER ::  &
           blas_upper_triangular = blas_symmetry_type(6)
      TYPE (blas_field_type), PARAMETER :: blas_complex = blas_field_type(1)
      TYPE (blas_field_type), PARAMETER :: blas_real = blas_field_type(2)
      TYPE (blas_field_type), PARAMETER ::  &
           blas_double_precision = blas_field_type(3)
      TYPE (blas_field_type), PARAMETER ::  &
           blas_single_precision = blas_field_type(4)
      TYPE (blas_size_type), PARAMETER :: blas_num_rows = blas_field_type(1)
      TYPE (blas_size_type), PARAMETER :: blas_num_cols = blas_field_type(2)
      TYPE (blas_size_type), PARAMETER :: blas_num_nonzeros = blas_field_type(3)
      TYPE (blas_valid_handle_type), PARAMETER ::  &
           blas_valid_handle = blas_valid_handle_type(1)
      TYPE (blas_sparsity_optimization_type), PARAMETER ::  &
           blas_regular = blas_sparsity_optimization_type(1)
      TYPE (blas_sparsity_optimization_type), PARAMETER ::  &
           blas_irregular = blas_sparsity_optimization_type(2)
      TYPE (blas_sparsity_optimization_type), PARAMETER ::  &
           blas_block = blas_sparsity_optimization_type(3)
      TYPE (blas_sparsity_optimization_type), PARAMETER ::  &
           blas_unassembled = blas_sparsity_optimization_type(4)
    END MODULE blas_operator_arguments
