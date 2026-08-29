MODULE AcousticSolver

USE Kinds, ONLY: WP, I32, SPEED_OF_SOUND, RHO_WATER
USE Kinds, ONLY: PREF => P_REF, IN_CLUSTER
IMPLICIT NONE

PRIVATE
PUBLIC :: AcousticSolver_t


! ------------------
! Abstract Base Type
! ------------------
TYPE, ABSTRACT :: AcousticSolver_t
    
    REAL(WP)          :: c       = SPEED_OF_SOUND       ! [m/s] Speed of sound in fluid
    REAL(WP)          :: rho     = RHO_WATER            ! [kg/m^3] Density of the fluid
    REAL(WP)          :: eps     = 1e-12_WP             ! [-] Small parameter for regularization
    REAL(WP)          :: p_ref   = PREF                 ! [Pa] Reference pressure
    LOGICAL           :: cluster = IN_CLUSTER           ! [-] Flag for cluster execution
    CHARACTER(len=50) :: solver_name = 'None'           ! [-] Solver name

    CONTAINS

    ! --- Deferred(abstract) procedures --- !
    PROCEDURE(get_name_i)        , DEFERRED :: get_name
    PROCEDURE(compute_pressure_i), DEFERRED :: compute_pressure

END TYPE AcousticSolver_t


! --------------------------------------------------
! Abstract interfaces for deferred(abstract) methods
! --------------------------------------------------
! Get solver name
ABSTRACT INTERFACE
    FUNCTION get_name_i(self) RESULT(name)
        IMPORT :: AcousticSolver_t
        CLASS(AcousticSolver_t), INTENT(INOUT) :: self
        CHARACTER(len=:), ALLOCATABLE          :: name
    END FUNCTION get_name_i
END INTERFACE


! Compute pressure
ABSTRACT INTERFACE
    SUBROUTINE compute_pressure_i(self, observers, block_size, total_pressure)
        IMPORT AcousticSolver_t, WP, I32
        CLASS(AcousticSolver_t) , INTENT(INOUT)        :: self
        REAL(WP)                , INTENT(IN)           :: observers(:,:)        ! [m] Observers coords (Nobs,3)
        INTEGER(I32)            , INTENT(IN), OPTIONAL :: block_size            ! [Pa] Pressure field (Nfreqs, Nobs)
        COMPLEX(WP), ALLOCATABLE, INTENT(OUT)          :: total_pressure(:,:)   ! [-] Number of observers to process simultaneously
    END SUBROUTINE compute_pressure_i
END INTERFACE


END MODULE AcousticSolver