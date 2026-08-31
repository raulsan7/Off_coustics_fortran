MODULE Kinds

IMPLICIT NONE
PRIVATE
PUBLIC :: I32, I64, WP, PI, I1
PUBLIC :: RHO_WATER, SPEED_OF_SOUND, P_REF
PUBLIC :: IN_CLUSTER

! ----------------------------------------
! Integer precision kinds
! ----------------------------------------
INTEGER, PARAMETER :: I32 = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: I64 = SELECTED_INT_KIND(18)


! ----------------------------------------
! Real precision kinds
! ----------------------------------------
INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6, 37)
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15, 307)
INTEGER, PARAMETER :: WP = DP                           ! Working precision, change when desired


! ----------------------------------------
! Physical and mathematical constants
! ----------------------------------------
REAL(WP)   , PARAMETER :: PI = 3.14159265358979323846_WP
COMPLEX(WP), PARAMETER :: I1 = (0.0_WP, 1.0_WP)         ! Imaginary unit

! Acoustic default properties of seawater
REAL(WP), PARAMETER :: RHO_WATER      = 1025.0_WP       ! [kg/m^3] Density
REAL(WP), PARAMETER :: SPEED_OF_SOUND = 1500.0_WP       ! [m/s] Sound velocity
REAL(WP), PARAMETER :: P_REF          = 1.0E-6_WP       ! [Pa] Reference pressure

! ----------------------------------------
! Other parameters
! ----------------------------------------
LOGICAL , PARAMETER :: IN_CLUSTER = .true.             ! [-] Wheter we are runnning in a cluster or not


END MODULE Kinds