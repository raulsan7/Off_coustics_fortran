PROGRAM MAIN

USE Kinds
USE IOUtils
USE WindTurbine
USE TurbineTypes

IMPLICIT NONE

! Local variables
CHARACTER(len=*), PARAMETER :: SD_path = "../OP_output/DTU_DeltaWind_mn_ws11.4.SD.sum.yaml"
CHARACTER(len=*), PARAMETER :: OP_path = "../OP_output/DTU_DeltaWind_mn_ws11.4.outb"
INTEGER(I32)    , PARAMETER :: Nm = 8, Nn=5

TYPE(DTU10MWMonopile) :: monopile

REAL(WP) :: start_time, end_time, elapsed_time
REAL(WP) :: elapsed_display
CHARACTER(len=8) :: tag


CALL cpu_time(start_time)
! ---------------------------------------------------------
! 1. Initialization and Banner
! ---------------------------------------------------------
CALL monopile % init( &
        rootname  = "DTU_DeltaWind_mn_ws11.4", &
        output_dir= "../OP_output/", &
        save_dir  = "./turbine_acoustic_data/", &
        WindSpeed = 11.4_dp, &
        WindDir   = 45.0_dp, &
        Depth     = 30.0_dp, &
        Nmembers  = 8, &
        Nnodes    = 5, &
        debug    = .true. )

CALL monopile % read_input(verbose=.true.)

! print *, "==================================================="
! print *, " OFF-Coustics: Offshore Acoustic Simulator         "
! print *, " Version: 0.1.0 (HPC Fortran Edition)              "
! print *, "==================================================="
! print *, "-> Precision (wp) initialized."
! print *, "-> Environment:"
! print '(A, F0.2, A)', "   * Water Density : ", RHO_WATER, " kg/m^3"
! print '(A, F0.2, A)', "   * Sound Speed   : ", SPEED_OF_SOUND, " m/s"
! print *, "---------------------------------------------------"

! ---------------------------------------------------------
! 2. Geometry and Kinematics Loading
! ---------------------------------------------------------


! ---------------------------------------------------------
! 3. Acoustic Solver Execution (To be implemented)
! ---------------------------------------------------------
! print *, "-> Assembling acoustic solvers..."
! TODO: initialize method_images or analytical_normal_modes
! TODO: call solver%compute_pressure_field()


! ---------------------------------------------------------
! 4. Results Output and Cleanup (To be implemented)
! ---------------------------------------------------------
! print *, "-> Exporting HDF5 results..."
! TODO: call h5_exporter%write_results()

CALL format_elapsed(start_time, elapsed_display, tag)

print*,' '
print *, "==================================================="
print '(A, F8.3, A)', '    Total elapsed time = ', elapsed_display, ' (' // trim(tag) // ')'
print *, "==================================================="

END PROGRAM MAIN