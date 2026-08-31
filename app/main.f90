! RUN COMMAND: 
! fpm clean                     
! fpm run --flag "-O3 -fopenmp" (only if IN_CLUSTER = .true. in Kinds.f90 else fpm run)

PROGRAM MAIN

USE Kinds, ONLY: I32, WP
USE TurbineTypes, ONLY: DTU10MWMonopile
USE MethodImages, ONLY: MethodImages_t
USE MathUtils, ONLY: format_elapsed

IMPLICIT NONE

! Local variables
CHARACTER(len=*), PARAMETER :: SD_path = "../OP_output/DTU_DeltaWind_mn_ws11.4.SD.sum.yaml"
CHARACTER(len=*), PARAMETER :: OP_path = "../OP_output/DTU_DeltaWind_mn_ws11.4.outb"
INTEGER(I32)    , PARAMETER :: Nm = 8, Nn = 5, unit_num = 10

TYPE(DTU10MWMonopile) :: monopile
TYPE(MethodImages_t)  :: acoustic_model

REAL(WP) :: start_time, elapsed_display
CHARACTER(len=8) :: tag


CALL cpu_time(start_time)
! ---------------------------------------------------------
! 1. Initialization and Banner
! ---------------------------------------------------------
CALL monopile % init( &
        rootname  = "DTU_DeltaWind_mn_ws11.4", &
        output_dir= "../OP_output/", &
        save_dir  = "./turbine_acoustic_data/", &
        WindSpeed = 11.4_WP, &
        WindDir   = 0.0_WP, &
        Depth     = 30.0_WP, &
        Nmembers  = 8, &
        Nnodes    = 5, &
        debug    = .true. )

CALL monopile % read_input(verbose=.true.)
CALL monopile % compute_force(filter_freqs=.true., verbose=.true.)
CALL acoustic_model % init(monopile, debug = .true.)

CALL acoustic_model % run_all()


! ---------- DISPLAY ELPASED TIME ---------- !
CALL format_elapsed(start_time, elapsed_display, tag)

print*,' '
print *, "==================================================="
print '(A, F8.3, A)', '    Total elapsed time = ', elapsed_display, ' (' // trim(tag) // ')'
print *, "==================================================="

END PROGRAM MAIN