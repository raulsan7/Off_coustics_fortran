! RUN COMMAND: 
! fpm clean                     
! fpm run --flag "-O3 -fopenmp" (only if IN_CLUSTER = .true. in Kinds.f90 else fpm run)

PROGRAM MAIN

USE Kinds, ONLY: I32, WP
USE MathUtils, ONLY: format_elapsed
USE MethodImages, ONLY: MethodImages_t
USE TurbineTypes, ONLY: DTU10MWMonopile, DTU10MWFloating
USE AnalyticalNormalModes, ONLY: AnalyticalNormalModes_t

IMPLICIT NONE

! Local variables
CHARACTER(len=*), PARAMETER :: SD_path = "../OP_output/DTU_DeltaWind_mn_ws11.4.SD.sum.yaml"
CHARACTER(len=*), PARAMETER :: OP_path = "../OP_output/DTU_DeltaWind_mn_ws11.4.outb"
INTEGER(I32)    , PARAMETER :: Nm = 8, Nn = 5, unit_num = 10

TYPE(DTU10MWMonopile) :: monopile
TYPE(DTU10MWFloating) :: floating
TYPE(MethodImages_t)  :: acoustic_model
TYPE(AnalyticalNormalModes_t) :: analytical_model

REAL(WP) :: start_time, elapsed_display
CHARACTER(len=8) :: tag


CALL cpu_time(start_time)
! ---------------------------------------------------------
! 1. Initialization and Banner
! ---------------------------------------------------------
! CALL monopile % init( &
!         rootname  = "DTU_DeltaWind_mn_ws11.4", &
!         output_dir= "../OP_output/", &
!         save_dir  = "./turbine_acoustic_data/", &
!         WindSpeed = 11.4_WP, &
!         WindDir   = 0.0_WP, &
!         Depth     = 30.0_WP, &
!         Nmembers  = 8, &
!         Nnodes    = 5)

! CALL monopile % read_input(verbose=.true.)
! CALL monopile % compute_force(filter_freqs=.true., verbose=.true.)
! CALL acoustic_model % init(monopile, debug = .true., name="plot_mn_SD30")
! CALL acoustic_model % run_all()
! CALL analytical_model % init(monopile, verbose=.true., name="plot_mn_ANM", debug=.true.)
! CALL analytical_model % run_all()
! CALL analytical_model % run_polar()

CALL floating % init( &
        rootname  = "DTU_DeltaWind_fl_ws11.4", &
        output_dir= "../OP_output/", &
        save_dir  = "./turbine_acoustic_data/", &
        WindSpeed = 11.4_WP, &
        WindDir   = 0.0_WP, &
        Depth     = 350.0_WP, &
        Nmembers  = 9, &
        Nnodes    = 9)

CALL floating % read_input(verbose=.true.)
CALL floating % compute_force(filter_freqs=.true., verbose=.true.)
CALL acoustic_model % init(floating, debug = .true., name="plot_fl_SD30", Lower_HBC=-350.0_WP)
! CALL acoustic_model % run_all()
CALL acoustic_model % run_spectrums(z_obs = 15.0_WP)
CALL acoustic_model % run_polar(z = -15.0_WP)
CALL acoustic_model % run_decay(z=-15.0_WP)
CALL acoustic_model % run_line()
CALL acoustic_model % run_sphere()
CALL acoustic_model % run_sliceXY(z = -15.0_WP)
CALL acoustic_model % run_sliceXZ()
CALL acoustic_model % run_cylinder()


! ---------- DISPLAY ELPASED TIME ---------- !
CALL format_elapsed(start_time, elapsed_display, tag)

print*,' '
print *, "==================================================="
print '(A, F8.3, A)', '    Total elapsed time = ', elapsed_display, ' (' // trim(tag) // ')'
print *, "==================================================="

END PROGRAM MAIN