MODULE TurbineTypes

USE Kinds
USE WindTurbine

IMPLICIT NONE 
PRIVATE

PUBLIC:: DTU10MWMonopile

TYPE, EXTENDS(WindTurbine_t) :: DTU10MWMonopile
    ! Physical and structural parameters specific to the DTU 10 MW Monopile
    REAL(WP) :: D        = 9.0_dp      ! [m] Base monopile outer diameter
    REAL(WP) :: rho_wat  = 1025.0_dp   ! [kg/m^3] Seawater density
    REAL(WP) :: rho_mat  = 8500.0_dp   ! [kg/m^3] Material density (steel)
    REAL(WP) :: wet_area = 848.12_dp   ! [m^2] Wetted area

    ! Path to rotor speed curve CSV file
    CHARACTER(len=:), ALLOCATABLE :: path_rpm

    CONTAINS
    PROCEDURE :: init                          => init_monopile
    ! PROCEDURE :: compute_force                 => compute_force_monopile
    ! PROCEDURE :: compute_masss                 => compute_mass_monopile
    ! PROCEDURE :: filter_frequencies            => filter_frequencies_monopile
    ! PROCEDURE :: get_impedance_corrected_force => get_impedance_corrected_force_monopile


END TYPE DTU10MWMonopile


CONTAINS

SUBROUTINE init_monopile(self, debug, rootname, output_dir, save_dir, save_name, &
                         WindSpeed, WindDir, Depth, AxisPos, BariPos, Binary, &
                         Nmembers, Nnodes)
    CLASS(DTU10MWMonopile), INTENT(INOUT) :: self
    LOGICAL               , INTENT(IN), OPTIONAL :: debug
    CHARACTER(len=*)      , INTENT(IN), OPTIONAL :: rootname
    CHARACTER(len=*)      , INTENT(IN), OPTIONAL :: output_dir
    CHARACTER(len=*)      , INTENT(IN), OPTIONAL :: save_dir
    CHARACTER(len=*)      , INTENT(IN), OPTIONAL :: save_name
    REAL(WP)              , INTENT(IN), OPTIONAL :: WindSpeed
    REAL(WP)              , INTENT(IN), OPTIONAL :: WindDir
    REAL(WP)              , INTENT(IN), OPTIONAL :: Depth
    REAL(WP)              , INTENT(IN), OPTIONAL :: AxisPos(2)
    REAL(WP)              , INTENT(IN), OPTIONAL :: BariPos(2)
    LOGICAL               , INTENT(IN), OPTIONAL :: Binary
    INTEGER(I32)          , INTENT(IN), OPTIONAL :: Nmembers
    INTEGER(I32)          , INTENT(IN), OPTIONAL :: Nnodes

    ! Local variables
    LOGICAL :: exists = .false.

    CALL init( self, &
            debug      = debug,      &
            rootname   = rootname,   &
            output_dir = output_dir, &
            save_dir   = save_dir,   &
            save_name  = save_name,  &
            WindSpeed  = WindSpeed,  &
            WindDir    = WindDir,    &
            Depth      = Depth,      &
            AxisPos    = AxisPos,    &
            BariPos    = BariPos,    &
            Binary     = Binary,     &
            Nmembers   = Nmembers,   &
            Nnodes     = Nnodes      )

    ! Set monopile-specific physical parameters
    self % D        = 9.0_WP
    self % rho_wat  = 1025.0_WP
    self % rho_mat  = 8500._WP
    self % wet_area = 848.12_WP

    ! Define path to rpm curve
    self % path_rpm = "../wind_speed_curves_DTU_10MW/rpm_ws.csv"

    ! Verify that the rpm curve file exists
    inquire(file= self % path_rpm, exist=exists)
    if (.not. exists) error stop "DTU10MWMonopile % init: RPM curve file not found: " // self % path_rpm

    ! Set the case type
    self % case_type = "Monopile"

END SUBROUTINE init_monopile



END MODULE TurbineTypes