PROGRAM MAIN

USE Kinds
USE IOUtils

IMPLICIT NONE

! Local variables
CHARACTER(len=*), PARAMETER :: SD_path = "../OP_output/DTU_DeltaWind_mn_ws11.4.SD.sum.yaml"
CHARACTER(len=*), PARAMETER :: OP_path = "../OP_output/DTU_DeltaWind_mn_ws11.4.outb"
INTEGER(I32), PARAMETER :: Nm = 8, Nn=5

! Variables for SubDyn node geometry
REAL(WP), ALLOCATABLE :: Nodes(:,:,:)

! Variables to receive READ_INPUT_SD outputs
REAL(WP), ALLOCATABLE :: time_array(:)
REAL(WP), ALLOCATABLE :: sd_tensor(:,:,:,:)
CHARACTER(len=32)     :: unit_label
INTEGER(I32)          :: status

! INTEGER(I32) :: i, j

! ---------------------------------------------------------
! 1. Initialization and Banner
! ---------------------------------------------------------
print *, "==================================================="
print *, " OFF-Coustics: Offshore Acoustic Simulator         "
print *, " Version: 0.1.0 (HPC Fortran Edition)              "
print *, "==================================================="
print *, "-> Precision (wp) initialized."
print *, "-> Environment:"
print '(A, F0.2, A)', "   * Water Density : ", RHO_WATER, " kg/m^3"
print '(A, F0.2, A)', "   * Sound Speed   : ", SPEED_OF_SOUND, " m/s"
print *, "---------------------------------------------------"

! ---------------------------------------------------------
! 2. Geometry and Kinematics Loading
! ---------------------------------------------------------
print *, "-> Parsing SubDyn summary..."
ALLOCATE(Nodes(Nm,Nn,3), source=0.0_WP)
CALL GET_SDSUM_VARIABLES(SD_path, Nmembers=Nm, Nnodes=Nn, verbose=.true., Nodes=Nodes)

print *, "-> Reading SubDyn kinematics from OpenFAST binary..."
CALL READ_INPUT_SD(filename=OP_path, Nmembers=Nm, Nnodes=Nn, &
                   verbose=.true., Time_out=time_array, Array_out=sd_tensor,      &
                   unit_out=unit_label, status=status)




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

print *, "==================================================="
print *, " Simulation completed successfully.                "
print *, "==================================================="

END PROGRAM MAIN