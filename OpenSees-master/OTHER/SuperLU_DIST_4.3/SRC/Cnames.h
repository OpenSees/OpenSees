/*! @file
 * \brief Macro definitions
 *
 * <pre>
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 * </pre>
 */

#ifndef __SUPERLU_CNAMES /* allow multiple inclusions */
#define __SUPERLU_CNAMES

/*
 * These macros define how C routines will be called.  ADD_ assumes that
 * they will be called by fortran, which expects C routines to have an
 * underscore postfixed to the name (Suns, and the Intel expect this).
 * NOCHANGE indicates that fortran will be calling, and that it expects
 * the name called by fortran to be identical to that compiled by the C
 * (RS6K's do this).  UPCASE says it expects C routines called by fortran
 * to be in all upcase (CRAY wants this). 
 */

#define ADD_       0
#define NOCHANGE   1
#define UPCASE     2
#define C_CALL     3

#ifdef UpCase
#define F77_CALL_C UPCASE
#endif

#ifdef NoChange
#define F77_CALL_C NOCHANGE
#endif

#ifdef Add_
#define F77_CALL_C ADD_
#endif

#ifndef F77_CALL_C
#define F77_CALL_C ADD_
#endif

#if (F77_CALL_C == ADD_)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine
 * No redefinition necessary to have following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void dgemm_(...)
 *
 * This is the default.
 */
/* These are the functions defined in F90 wraper */
#define f_create_gridinfo_handle         f_create_gridinfo_handle_
#define f_create_options_handle          f_create_options_handle_
#define f_create_ScalePerm_handle        f_create_scaleperm_handle_
#define f_create_LUstruct_handle         f_create_lustruct_handle_
#define f_create_SOLVEstruct_handle      f_create_solvestruct_handle_
#define f_create_SuperMatrix_handle      f_create_supermatrix_handle_
#define f_destroy_gridinfo_handle        f_destroy_gridinfo_handle_
#define f_destroy_options_handle         f_destroy_options_handle_
#define f_destroy_ScalePerm_handle       f_destroy_scaleperm_handle_
#define f_destroy_LUstruct_handle        f_destroy_lustruct_handle_
#define f_destroy_SOLVEstruct_handle     f_destroy_solvestruct_handle_
#define f_destroy_SuperMatrix_handle     f_destroy_supermatrix_handle_
#define f_create_SuperLUStat_handle      f_create_superlustat_handle_
#define f_destroy_SuperLUStat_handle     f_destroy_superlustat_handle_
#define f_get_gridinfo                   f_get_gridinfo_
#define f_get_SuperMatrix                f_get_supermatrix_
#define f_set_SuperMatrix                f_set_supermatrix_
#define f_get_CompRowLoc_Matrix          f_get_comprowloc_matrix_ 
#define f_set_CompRowLoc_Matrix          f_set_comprowloc_matrix_
#define f_get_superlu_options            f_get_superlu_options_
#define f_set_superlu_options            f_set_superlu_options_
#define f_set_default_options            f_set_default_options_
#define f_superlu_gridinit               f_superlu_gridinit_
#define f_superlu_gridmap                f_superlu_gridmap_
#define f_superlu_gridexit               f_superlu_gridexit_
#define f_ScalePermstructInit            f_scalepermstructinit_
#define f_ScalePermstructFree            f_scalepermstructfree_
#define f_PStatInit                      f_pstatinit_
#define f_PStatFree                      f_pstatfree_
#define f_LUstructInit                   f_lustructinit_
#define f_LUstructFree                   f_lustructfree_
#define f_Destroy_LU                     f_destroy_lu_
#define f_dCreate_CompRowLoc_Mat_dist    f_dcreate_comprowloc_mat_dist_
#define f_zCreate_CompRowLoc_Mat_dist    f_zcreate_comprowloc_mat_dist_
#define f_Destroy_CompRowLoc_Mat_dist    f_destroy_comprowloc_mat_dist_
#define f_Destroy_SuperMat_Store_dist    f_destroy_supermat_store_dist_
#define f_dSolveFinalize                 f_dsolvefinalize_
#define f_zSolveFinalize                 f_zsolvefinalize_
#define f_pdgssvx                        f_pdgssvx_
#define f_pzgssvx                        f_pzgssvx_
#define f_dcreate_dist_matrix            f_dcreate_dist_matrix_
#define f_zcreate_dist_matrix            f_zcreate_dist_matrix_
#define f_check_malloc                   f_check_malloc_
#endif

#if (F77_CALL_C == UPCASE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void DGEMM(...)
 */
/* BLAS */
#define sasum_    SASUM
#define isamax_   ISAMAX
#define scopy_    SCOPY
#define sscal_    SSCAL
#define sger_     SGER
#define snrm2_    SNRM2
#define ssymv_    SSYMV
#define sdot_     SDOT
#define saxpy_    SAXPY
#define ssyr2_    SSYR2
#define srot_     SROT
#define sgemv_    SGEMV
#define strsv_    STRSV
#define sgemm_    SGEMM
#define strsm_    STRSM

#define dasum_    DASUM
#define idamax_   IDAMAX
#define dcopy_    DCOPY
#define dscal_    DSCAL
#define dger_     DGER
#define dnrm2_    DNRM2
#define dsymv_    DSYMV
#define ddot_     DDOT
#define daxpy_    DAXPY
#define dsyr2_    DSYR2
#define drot_     DROT
#define dgemv_    DGEMV
#define dtrsv_    DTRSV
#define dgemm_    DGEMM
#define dtrsm_    DTRSM

#define scasum_   SCASUM
#define icamax_   ICAMAX
#define ccopy_    CCOPY
#define cscal_    CSCAL
#define scnrm2_   SCNRM2
#define caxpy_    CAXPY
#define cgemv_    CGEMV
#define ctrsv_    CTRSV
#define cgemm_    CGEMM
#define ctrsm_    CTRSM
#define cgerc_    CGERC
#define chemv_    CHEMV
#define cher2_    CHER2

#define dzasum_   DZASUM
#define izamax_   IZAMAX
#define zcopy_    ZCOPY
#define zscal_    ZSCAL
#define dznrm2_   DZNRM2
#define zaxpy_    ZAXPY
#define zgemv_    ZGEMV
#define ztrsv_    ZTRSV
#define zgemm_    ZGEMM
#define ztrsm_    ZTRSM
#define zgerc_    ZGERC
#define zhemv_    ZHEMV
#define zher2_    ZHER2
#define zgeru_    ZGERU

/*
#define mc64id_dist     MC64ID_DIST
#define mc64ad_dist     MC64AD_DIST
*/
#define c_bridge_dgssv_               C_BRIDGE_DGSSV
#define c_fortran_slugrid_            C_FORTRAN_SLUGRID
#define c_fortran_pdgssvx_            C_FORTRAN_PDGSSVX
#define c_fortran_pdgssvx_ABglobal_   C_FORTRAN_PDGSSVX_ABGLOBAL
#define c_fortran_pzgssvx_            C_FORTRAN_PZGSSVX
#define c_fortran_pzgssvx_ABglobal_   C_FORTRAN_PZGSSVX_ABGLOBAL

/* These are the functions defined in F90 wraper */
#define f_create_gridinfo_handle         F_CREATE_GRIDINFO_HANDLE
#define f_create_options_handle          F_CREATE_OPTIONS_HANDLE
#define f_create_ScalePerm_handle        F_CREATE_SCALEPERM_HANDLE
#define f_create_LUstruct_handle         F_CREATE_LUSTRUCT_HANDLE
#define f_create_SOLVEstruct_handle      F_CREATE_SOLVESTRUCT_HANDLE
#define f_create_SuperMatrix_handle      F_CREATE_SUPERMATRIX_HANDLE
#define f_destroy_gridinfo_handle        F_DESTROY_GRIDINFO_HANDLE
#define f_destroy_options_handle         F_DESTROY_OPTIONS_HANDLE
#define f_destroy_ScalePerm_handle       F_DESTROY_SCALEPERM_HANDLE
#define f_destroy_LUstruct_handle        F_DESTROY_LUSTRUCT_HANDLE
#define f_destroy_SOLVEstruct_handle     F_DESTROY_SOLVESTRUCT_HANDLE
#define f_destroy_SuperMatrix_handle     F_DESTROY_SUPERMATRIX_HANDLE
#define f_create_SuperLUStat_handle      F_CREATE_SUPERLUSTAT_HANDLE
#define f_destroy_SuperLUStat_handle     F_DESTROY_SUPERLUSTAT_HANDLE
#define f_get_gridinfo                   F_GET_GRIDINFO
#define f_get_SuperMatrix                F_GET_SUPERMATRIX
#define f_set_SuperMatrix                F_SET_SUPERMATRIX
#define f_get_CompRowLoc_Matrix          F_GET_COMPROWLOC_MATRIX
#define f_set_CompRowLoc_Matrix          F_SET_COMPROWLOC_MATRIX
#define f_get_superlu_options            F_GET_SUPERLU_OPTIONS
#define f_set_superlu_options            F_SET_SUPERLU_OPTIONS
#define f_set_default_options            F_SET_DEFAULT_OPTIONS
#define f_superlu_gridinit               F_SUPERLU_GRIDINIT
#define f_superlu_gridmap                F_SUPERLU_GRIDMAP
#define f_superlu_gridexit               F_SUPERLU_GRIDEXIT
#define f_ScalePermstructInit            F_SCALEPERMSTRUCTINIT
#define f_ScalePermstructFree            F_SCALEPERMSTRUCTFREE
#define f_PStatInit                      F_PSTATINIT
#define f_PStatFree                      F_PSTATFREE
#define f_LUstructInit                   F_LUSTRUCTINIT
#define f_LUstructFree                   F_LUSTRUCTFREE
#define f_Destroy_LU                     F_DESTROY_LU
#define f_dCreate_CompRowLoc_Mat_dist    F_DCREATE_COMPROWLOC_MAT_DIST
#define f_zCreate_CompRowLoc_Mat_dist    F_ZCREATE_COMPROWLOC_MAT_DIST
#define f_Destroy_CompRowLoc_Mat_dist    F_DESTROY_COMPROWLOC_MAT_DIST
#define f_Destroy_SuperMat_Store_dist    F_DESTROY_SUPERMAT_STORE_DIST
#define f_dSolveFinalize                 F_DSOLVEFINALIZE
#define f_zSolveFinalize                 F_ZSOLVEFINALIZE
#define f_pdgssvx                        F_PDGSSVX
#define f_pzgssvx                        F_PZGSSVX
#define f_dcreate_dist_matrix            F_DCREATE_DIST_MATRIX
#define f_zcreate_dist_matrix            F_ZCREATE_DIST_MATRIX
#define f_check_malloc                   F_CHECK_MALLOC
#endif

#if (F77_CALL_C == NOCHANGE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * for following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void dgemm(...)
 */
/* BLAS */
#define sasum_    sasum
#define isamax_   isamax
#define scopy_    scopy
#define sscal_    sscal
#define sger_     sger
#define snrm2_    snrm2
#define ssymv_    ssymv
#define sdot_     sdot
#define saxpy_    saxpy
#define ssyr2_    ssyr2
#define srot_     srot
#define sgemv_    sgemv
#define strsv_    strsv
#define sgemm_    sgemm
#define strsm_    strsm

#define dasum_    dasum
#define idamax_   idamax
#define dcopy_    dcopy
#define dscal_    dscal
#define dger_     dger
#define dnrm2_    dnrm2
#define dsymv_    dsymv
#define ddot_     ddot
#define daxpy_    daxpy
#define dsyr2_    dsyr2
#define drot_     drot
#define dgemv_    dgemv
#define dtrsv_    dtrsv
#define dgemm_    dgemm
#define dtrsm_    dtrsm

#define scasum_   scasum
#define icamax_   icamax
#define ccopy_    ccopy
#define cscal_    cscal
#define scnrm2_   scnrm2
#define caxpy_    caxpy
#define cgemv_    cgemv
#define ctrsv_    ctrsv
#define cgemm_    cgemm
#define ctrsm_    ctrsm
#define cgerc_    cgerc
#define chemv_    chemv
#define cher2_    cher2

#define dzasum_   dzasum
#define izamax_   izamax
#define zcopy_    zcopy
#define zscal_    zscal
#define dznrm2_   dznrm2
#define zaxpy_    zaxpy
#define zgemv_    zgemv
#define ztrsv_    ztrsv
#define zgemm_    zgemm
#define ztrsm_    ztrsm
#define zgerc_    zgerc
#define zhemv_    zhemv
#define zher2_    zher2
#define zgeru_    zgeru

/*
#define mc64id_dist         mc64id_dist
#define mc64ad_dist         mc64ad_dist
*/

#define c_bridge_dgssv_               c_bridge_dgssv
#define c_fortran_slugrid_            c_fortran_slugrid
#define c_fortran_pdgssvx_            c_fortran_pdgssvx
#define c_fortran_pdgssvx_ABglobal_   c_fortran_pdgssvx_abglobal
#define c_fortran_pzgssvx_            c_fortran_pzgssvx
#define c_fortran_pzgssvx_ABglobal_   c_fortran_pzgssvx_abglobal

/* These are the functions defined in F90 wraper */
#define f_create_gridinfo_handle         f_create_gridinfo_handle
#define f_create_options_handle          f_create_options_handle
#define f_create_ScalePerm_handle        f_create_scaleperm_handle
#define f_create_LUstruct_handle         f_create_lustruct_handle
#define f_create_SOLVEstruct_handle      f_create_solvestruct_handle
#define f_create_SuperMatrix_handle      f_create_supermatrix_handle
#define f_destroy_gridinfo_handle        f_destroy_gridinfo_handle
#define f_destroy_options_handle         f_destroy_options_handle
#define f_destroy_ScalePerm_handle       f_destroy_scaleperm_handle
#define f_destroy_LUstruct_handle        f_destroy_lustruct_handle
#define f_destroy_SOLVEstruct_handle     f_destroy_solvestruct_handle
#define f_destroy_SuperMatrix_handle     f_destroy_supermatrix_handle
#define f_create_SuperLUStat_handle      f_create_superlustat_handle
#define f_destroy_SuperLUStat_handle     f_destroy_superlustat_handle
#define f_get_gridinfo                   f_get_gridinfo
#define f_get_SuperMatrix                f_get_supermatrix
#define f_set_SuperMatrix                f_set_supermatrix
#define f_get_CompRowLoc_Matrix          f_get_comprowloc_matrix 
#define f_set_CompRowLoc_Matrix          f_set_comprowloc_matrix
#define f_get_superlu_options            f_get_superlu_options
#define f_set_superlu_options            f_set_superlu_options
#define f_set_default_options            f_set_default_options
#define f_superlu_gridinit               f_superlu_gridinit
#define f_superlu_gridmap                f_superlu_gridmap
#define f_superlu_gridexit               f_superlu_gridexit
#define f_ScalePermstructInit            f_scalepermstructinit
#define f_ScalePermstructFree            f_scalepermstructfree
#define f_PStatInit                      f_pstatinit
#define f_PStatFree                      f_pstatfree
#define f_LUstructInit                   f_lustructinit
#define f_LUstructFree                   f_lustructfree
#define f_Destroy_LU                     f_destroy_lu
#define f_dCreate_CompRowLoc_Mat_dist    f_dcreate_comprowloc_mat_dist
#define f_Destroy_CompRowLoc_Mat_dist    f_destroy_comprowloc_mat_dist
#define f_Destroy_SuperMat_Store_dist    f_destroy_supermat_store_dist
#define f_dSolveFinalize                 f_dsolvefinalize
#define f_zSolveFinalize                 f_zsolvefinalize
#define f_pdgssvx                        f_pdgssvx
#define f_pzgssvx                        f_pzgssvx
#define f_dcreate_dist_matrix            f_dcreate_dist_matrix
#define f_zcreate_dist_matrix            f_zcreate_dist_matrix
#define f_check_malloc                   f_check_malloc
#endif

#endif /* __SUPERLU_CNAMES */
