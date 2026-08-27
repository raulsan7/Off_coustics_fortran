MODULE FFTW3

USE, INTRINSIC :: iso_c_binding
IMPLICIT NONE

PUBLIC

! FFTW Constants
INTEGER(C_INT), PARAMETER, PUBLIC :: FFTW_ESTIMATE   = 64
INTEGER(C_INT), PARAMETER, PUBLIC :: FFTW_MEASURE    = 0
INTEGER(C_INT), PARAMETER, PUBLIC :: FFTW_PATIENT    = 32
INTEGER(C_INT), PARAMETER, PUBLIC :: FFTW_EXHAUSTIVE = 8

! FFTW functions
INTERFACE
    TYPE(C_PTR) FUNCTION fftw_plan_many_dft_r2c(rank, n, howmany, &
                                                in, inembed, istride, idist, &
                                                out, onembed, ostride, odist, &
                                                flags) BIND(c, name='fftw_plan_many_dft_r2c')

        IMPORT :: C_INT, C_PTR, C_DOUBLE, C_DOUBLE_COMPLEX
        INTEGER(C_INT), VALUE :: rank
        INTEGER(C_INT), DIMENSION(*) :: n
        INTEGER(C_INT), VALUE :: howmany
        REAL(C_DOUBLE), DIMENSION(*) :: in
        INTEGER(C_INT), DIMENSION(*) :: inembed
        INTEGER(C_INT), VALUE :: istride, idist 
        COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(*) :: out
        INTEGER(C_INT), DIMENSION(*) :: onembed
        INTEGER(C_INT), VALUE :: ostride, odist
        INTEGER(C_INT), VALUE :: flags

    END FUNCTION fftw_plan_many_dft_r2c

    SUBROUTINE fftw_execute_dft_r2c(plan, in, out) BIND(c, name='fftw_execute_dft_r2c')

        IMPORT :: C_PTR, C_DOUBLE, c_double_complex
        TYPE(C_PTR), VALUE :: plan
        REAL(C_DOUBLE), DIMENSION(*) :: in
        COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(*):: out
    
    END SUBROUTINE fftw_execute_dft_r2c

    SUBROUTINE fftw_destroy_plan(plan) bind(c, name='fftw_destroy_plan')

        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: plan

    END SUBROUTINE fftw_destroy_plan

END INTERFACE


END MODULE FFTW3