MODULE WindTurbine

    
USE omp_lib
USE Kinds, ONLY: I32, WP, PI

IMPLICIT NONE
PRIVATE
PUBLIC :: WindTurbine_t, init

! ------------------
! Abstract base type
! ------------------
TYPE, ABSTRACT :: WindTurbine_t
    ! --- Input parameters --- !
    LOGICAL                       :: debug = .false.                ! [-] Using debug mode will print aditional data
    CHARACTER(len=:), ALLOCATABLE :: rootname                       ! [-] Name without extensions of the OpenFAST output files
    REAL(WP)                      :: WindSpeed = 0.0_WP             ! [m/s] Wind Speed in norm
    REAL(WP)                      :: WindDir = 0.0_WP               ! [deg] Wind direction 0 deg points to +x axis (anticlockwise from +x)
    REAL(WP)                      :: Depth = 0.0_WP                 ! [m] Water depth
    REAL(WP)                      :: AxisPos(2) = [0.0_WP, 0.0_WP]  ! [m] Position of the turbine axis in xy plane
    REAL(WP)                      :: BariPos(2) = [0.0_WP, 0.0_WP]  ! [m] Position of the turbine baricenter in xy plane
    INTEGER(I32)                  :: Nmembers = 0                   ! [-] Number of structural members in the OpenFAST model
    INTEGER(I32)                  :: Nnodes = 0                     ! [-] Number of structural nodes in the OpenFAST model

    ! --- Derived attirbutes --- !
    CHARACTER(len=:), ALLOCATABLE :: case_type                      ! [-] Tag that identifies case type
    LOGICAL                       :: in_farm = .false.              ! [-] Wheter current turbine is within a farm
    COMPLEX(WP), ALLOCATABLE      :: F(:,:,:)                       ! [N] Force array (Nfreqs, Nnodes_wet, 3)
    REAL(WP)   , ALLOCATABLE      :: Freqs(:)                       ! [Hz] Frequency array (Nfreqs)
    REAL(WP)   , ALLOCATABLE      :: x_all(:,:,:)                   ! [m] All nodes array (Nmembers, Nnodes, 3)
    REAL(WP)   , ALLOCATABLE      :: x(:,:)                         ! [m] Wetted nodes list (Nnodes_wet, 3)
    REAL(WP)   , ALLOCATABLE      :: time(:)                        ! [s] Time array
    REAL(WP)   , ALLOCATABLE      :: acc(:,:,:,:)                   ! [m/s^2] Acceleration array (NFreqs, Nmembers, Nnodes, 3)
    LOGICAL                       :: BariPos_set = .false.          ! [-] Wheter baricenter is inputed

    ! --- File paths ---
    CHARACTER(len=:), ALLOCATABLE :: outfile                        ! [-] OpenFast output file path .outb/.out
    CHARACTER(len=:), ALLOCATABLE :: SUM_SubDyn                     ! [-] SubDyn summary file path
    CHARACTER(len=:), ALLOCATABLE :: save_path                      ! [-] Acoustic output file


    CONTAINS
    ! --- Public type-bound procedures --- !
    PROCEDURE :: init                           ! Constructor
    PROCEDURE :: read_input                     ! Reads OpenFast outputs
    PROCEDURE :: translate_align_with_WindDir   ! Place nodes correctly


    ! --- Deferred(abstract) procedures --- ! Note: compute_force_i is only for the sign, compute_force is the name used in the subclases
    PROCEDURE(compute_force_i)                , DEFERRED :: compute_force
    PROCEDURE(compute_mass_i)                 , DEFERRED :: compute_mass
    PROCEDURE(filter_frequencies_i)           , DEFERRED :: filter_frequencies
    PROCEDURE(get_impedance_corrected_force_i), DEFERRED :: get_impedance_corrected_force

END TYPE WindTurbine_t


! --------------------------------------------------
! Abstract interfaces for deferred(abstract) methods
! --------------------------------------------------
! Compute force
ABSTRACT INTERFACE
    SUBROUTINE compute_force_i(self, filter_freqs, verbose, skipf)
        IMPORT :: WindTurbine_t, I32
        CLASS(WindTurbine_t), INTENT(INOUT) :: self
        LOGICAL             , INTENT(IN), OPTIONAL :: filter_freqs
        LOGICAL             , INTENT(IN), OPTIONAL :: verbose
        INTEGER(I32)        , INTENT(IN), OPTIONAL :: skipf
    END SUBROUTINE compute_force_i
END INTERFACE

! Compute mass
ABSTRACT INTERFACE
    SUBROUTINE compute_mass_i(self, mass, added_mass)
        IMPORT :: WindTurbine_t, WP
        CLASS(WindTurbine_t) , INTENT(INOUT):: self
        REAL(WP), ALLOCATABLE, INTENT(OUT) :: mass(:)
        REAL(WP), ALLOCATABLE, INTENT(OUT) :: added_mass(:)
    END SUBROUTINE compute_mass_i
END INTERFACE

! Filter frequencies
ABSTRACT INTERFACE
    FUNCTION filter_frequencies_i(self) RESULT(mask)
        IMPORT :: WindTurbine_t, WP
        CLASS(WindTurbine_t), INTENT(INOUT) :: self
        LOGICAL, ALLOCATABLE                :: mask(:)
    END FUNCTION filter_frequencies_i
END INTERFACE

! Impedance correction
ABSTRACT INTERFACE
    FUNCTION get_impedance_corrected_force_i(self, c) RESULT(corrected_force)
        IMPORT :: WindTurbine_t, WP
        CLASS(WindTurbine_t), INTENT(IN) :: self
        REAL(WP), INTENT(IN) :: c
        COMPLEX(WP), ALLOCATABLE :: corrected_force(:,:,:)
    END FUNCTION get_impedance_corrected_force_i
END INTERFACE


CONTAINS

SUBROUTINE init(self, debug, rootname, output_dir, save_dir, save_name, &
                WindSpeed, WindDir, Depth, AxisPos, BariPos, Binary, &
                Nmembers, Nnodes)
    CLASS(WindTurbine_t), INTENT(INOUT) :: self
    LOGICAL         , INTENT(IN), OPTIONAL :: debug
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: rootname
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: output_dir
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: save_dir
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: save_name
    REAL(WP)        , INTENT(IN), OPTIONAL :: WindSpeed
    REAL(WP)        , INTENT(IN), OPTIONAL :: WindDir
    REAL(WP)        , INTENT(IN), OPTIONAL :: Depth
    REAL(WP)        , INTENT(IN), OPTIONAL :: AxisPos(2)
    REAL(WP)        , INTENT(IN), OPTIONAL :: BariPos(2)
    LOGICAL         , INTENT(IN), OPTIONAL :: Binary
    INTEGER(I32)    , INTENT(IN), OPTIONAL :: Nmembers
    INTEGER(I32)    , INTENT(IN), OPTIONAL :: Nnodes

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: output_dir_, save_dir_, ext
    CHARACTER(len=:), ALLOCATABLE :: outfile_path, sum_path, save_path_temp


    ! Self assignment: scalar
    if (PRESENT(debug))     self%debug     = debug
    if (PRESENT(WindSpeed)) self%WindSpeed = WindSpeed
    if (PRESENT(WindDir))   self%WindDir   = WindDir
    if (PRESENT(Depth))     self%Depth     = Depth
    if (PRESENT(Nmembers))  self%Nmembers  = Nmembers
    if (PRESENT(Nnodes))    self%Nnodes    = Nnodes

    ! Self assignement: array
    if (PRESENT(AxisPos)) then
        self%AxisPos = AxisPos
    else
        self%AxisPos = [0.0_WP, 0.0_WP]
    end if

    if (PRESENT(BariPos)) then
        self%BariPos = BariPos
        self%Baripos_set = .true.
    else
        self%BariPos = [0.0_WP, 0.0_WP]
    end if

    ! Directory paths and file extension
    if (.not. PRESENT(rootname)) then
        error stop "WindTurbine_t $ init: rootname argument is required"
    end if
    self%rootname = rootname

    if (PRESENT(output_dir)) then
        output_dir_ = trim(output_dir)
    else
        output_dir_ = "../OP_output/"
    end if

    if (PRESENT(save_dir)) then
        save_dir_ = save_dir
    else 
        save_dir_ = "./turbine_acoustic_data/"
    end if

    if (PRESENT(Binary)) then
        if (Binary) then
            ext = ".outb"
        else 
            ext = ".out"
        end if
    else
        ext = ".outb"
    end if

    ! Build file paths
    outfile_path = trim(output_dir_) // trim(self%rootname) // trim(ext)
    sum_path     = trim(output_dir_) // trim(self%rootname) // ".SD.sum.yaml"

    if (PRESENT(save_name)) then
        save_path_temp = trim(save_dir_) // trim(save_name) // ".hdf5"
    else
        save_path_temp = trim(save_dir) ! If no save_name, store only the directory; save_acoustics will check suffix.
    end if

    self%outfile    = outfile_path
    self%SUM_SubDyn = sum_path
    self%save_path  = save_path_temp

    ! Other attibutes are left unallocated (case_type, Freqs, F ...)
    ! They will be set during read_input or compute_force

END SUBROUTINE init


SUBROUTINE read_input(self, in_farm, verbose, skip, From, Upto)
    USE IOUtils, ONLY: get_SDsum_variables, read_input_SD

    CLASS(WindTurbine_t), INTENT(INOUT)        :: self
    LOGICAL             , INTENT(IN), OPTIONAL :: in_farm
    LOGICAL             , INTENT(IN), OPTIONAL :: verbose
    INTEGER(I32)        , INTENT(IN), OPTIONAL :: skip
    REAL(WP)            , INTENT(IN), OPTIONAL :: From
    REAL(WP)            , INTENT(IN), OPTIONAL :: Upto

    ! Local variables
    LOGICAL :: in_farm_ = .false., verbose_ = .false., dir_exists
    INTEGER(I32) :: skip_ = 1, iostat, status
    REAL(WP) :: From_ = 0.0_WP, Upto_ = 1.0_WP, dt
    CHARACTER(len=:), ALLOCATABLE :: dir
    CHARACTER(len=16) :: unit

    ! Defaults
    if (PRESENT(in_farm)) in_farm_ = in_farm
    if (PRESENT(verbose)) verbose_ = verbose
    if (PRESENT(skip))    skip_    = skip
    if (PRESENT(From))    From_    = From
    if (PRESENT(Upto))    Upto_    = Upto
    if (PRESENT(in_farm)) in_farm_ = in_farm
    self % in_farm = in_farm_

    ! Create directories if needed
    if (len_trim(self % save_path) > 0) then
        ! Determine if save_path is a file or directory
        ! If it ends with '.hdf5' or '.npz', extract parent directory
        if (index(self % save_path, '.', back = .true.) > 0) then
            dir = self % save_path(1:index(self % save_path, '/', back = .true.)-1)
        else
            dir = trim(self % save_path)
        end if
        inquire(file=dir, exist=dir_exists)
        if (.not. dir_exists) then
            CALL execute_command_line('mkdir -p ' // trim(dir), exitstat=iostat)
            if (iostat /= 0)  error stop "WindTurbine % read_input: could not create directory " // trim(dir)
        end if
    end if

    ! Load initial positions from SubDyn summary file
    if (self % Nmembers > 0 .and. self % Nnodes >0) then
        CALL get_SDsum_variables(SD_path=self%SUM_SubDyn, Nmembers=self%Nmembers, &
                                 Nnodes=self%Nnodes, verbose=verbose_, Nodes=self%x_all)
        CALL read_input_SD(filename=self%outfile, what='acceleration', skip=skip_,Nmembers=self%Nmembers, &
                           Nnodes=self%Nnodes, From=From_, Upto=Upto_, verbose=verbose_, Time_out=self%time, &
                           Array_out=self%acc, unit_out=unit, status=status)
    else 
        error stop "WindTurbine % read_input: Nmembers and Nnodes must be positive integers to read SubDyn nodes"
    end if

    ! Compute baricenter if not provided
    if (.not. self % BariPos_set) then
        self % BariPos(1) = sum(self%x_all(:,:,1)) / real(self%Nmembers * self%Nnodes, WP)
        self % BariPos(2) = sum(self%x_all(:,:,2)) / real(self%Nmembers * self%Nnodes, WP)
        self % BariPos_set = .true.
    end if

    ! Align with wind direction and translate to axis position
    CALL self % translate_align_with_WindDir()

    ! Verbose summary
    if (verbose_) then
        dt = self%Time(2) - self%Time(1)
        print*,''
        print '(A)', repeat('=', 70)
        print '(A)', 'TURBINE: ' // trim(self%rootname)
        print '(A)', repeat('-', 68)
        print '(A, A)',  'Type: ', trim(self%case_type)
        print '(A, F5.1, A, F5.1, A)', 'Wind Speed: ', self%WindSpeed, ' m/s   | Wind Direction: ', self%WindDir, ' deg'
        print '(A, F4.1, A, F4.1, A, F4.1, A)', 'Water Depth: ', self%Depth,&
         ' m   | Position: (', self%AxisPos(1), ', ', self%AxisPos(2), ') m'
        print '(A, I6, A, F6.4, A)', 'Time Series: ', size(self%Time), ' samples @ ', dt, ' s/sample'
        print '(A)', repeat('=', 70)
    end if

    
END SUBROUTINE read_input


SUBROUTINE translate_align_with_WindDir(self)
    CLASS(WindTurbine_t), INTENT(INOUT) :: self

    REAL(WP) :: yaw, c, s
    REAL(WP) :: axis_xy(2)
    INTEGER(I32) :: i, j, m, n
    INTEGER(I32) :: nt, nm, nn
    REAL(WP) :: x_old, y_old

    yaw = self % WindDir * PI / 180.0_WP
    c   = cos(yaw)
    s   = sin(yaw)

    if (abs(self % WindDir) < 1.0e-12_WP) return

    axis_xy = self % AxisPos

    if (ALLOCATED(self % acc)) then
       nt = size(self % acc, 1)
       nm = size(self % acc, 2)
       nn = size(self % acc, 3)
       !$omp parallel do collapse(3)
       do n = 1, nn
          do m = 1, nm
             do i = 1, nt
                x_old = self % acc(i,m,n,1)
                y_old = self % acc(i,m,n,2)
                self % acc(i,m,n,1) = c * x_old - s * y_old
                self % acc(i,m,n,2) = s * x_old + c * y_old
             end do
          end do
       end do
       !$omp end parallel do
    end if

    if (ALLOCATED(self % x_all)) then
       nm = size(self % x_all, 1)
       nn = size(self % x_all, 2)
       !$omp parallel do collapse(2)
       do i = 1, nm
          do j = 1, nn
             x_old = self % x_all(i,j,1)
             y_old = self % x_all(i,j,2)
             self % x_all(i,j,1) = c * x_old - s * y_old + axis_xy(1)
             self % x_all(i,j,2) = s * x_old + c * y_old + axis_xy(2)
          end do
       end do
       !$omp end parallel do
    end if
END SUBROUTINE

END MODULE WindTurbine