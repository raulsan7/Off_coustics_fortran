MODULE WindTurbine

    
USE, INTRINSIC :: iso_fortran_env
USE Kinds, ONLY: I32, WP, PI
USE omp_lib
USE AcousticSolver

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

    ! --- Acoustic Solver placeholder --- !
    CLASS(AcousticSolver_t), ALLOCATABLE :: acoustic_solver

    ! --- File paths ---
    CHARACTER(len=:), ALLOCATABLE :: outfile                        ! [-] OpenFast output file path .outb/.out
    CHARACTER(len=:), ALLOCATABLE :: SUM_SubDyn                     ! [-] SubDyn summary file path
    CHARACTER(len=:), ALLOCATABLE :: save_path                      ! [-] Acoustic output file


    CONTAINS
    ! --- Public type-bound procedures --- !
    PROCEDURE :: init                           ! Constructor
    PROCEDURE :: read_input                     ! Reads OpenFast outputs
    PROCEDURE :: translate_align_with_WindDir   ! Place nodes correctly
    PROCEDURE :: set_acoustic_method            ! Selects which acoustic solver to use
    PROCEDURE :: save_parameters                ! Saves turbine parameters to save_path
    PROCEDURE :: check_acoustic_solver          ! Helper for run subroutines
    PROCEDURE :: check_observers_distances      ! Helper for run subroutines
    ! PROCEDURE :: save_acoustics                 ! Saves acoustics to save_path
    PROCEDURE :: run_spectrums                  ! Subroutines to compute acoustic pressure
    PROCEDURE :: run_polar
    PROCEDURE :: run_cylinder
    PROCEDURE :: run_decay
    PROCEDURE :: run_line
    PROCEDURE :: run_sliceXY
    PROCEDURE :: run_sliceXZ
    ! PROCEDURE :: run_sphere
    PROCEDURE :: run_all

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

SUBROUTINE set_acoustic_method(self, solver, unit_num)
    CLASS(WindTurbine_t)   , INTENT(INOUT) :: self
    CLASS(AcousticSolver_t), INTENT(IN) :: solver
    INTEGER(I32)           , INTENT(IN), OPTIONAL :: unit_num

    ! Local variables
    INTEGER(I32) :: unit_num_

    unit_num_ = 10; if (PRESENT(unit_num)) unit_num_ = unit_num

    if (ALLOCATED(self % acoustic_solver)) DEALLOCATE(self % acoustic_solver)
    ALLOCATE(self % acoustic_solver, source = solver)
    CALL save_parameters(self, unit_num_)

END SUBROUTINE set_acoustic_method


SUBROUTINE save_parameters(self, unit_num)
    CLASS(WindTurbine_t), INTENT(INOUT) :: self
    INTEGER(I32)        , INTENT(IN)    :: unit_num
    
    ! Save turbine parameters
    write(unit_num, *) "Rootname:", trim(self % rootname)
    write(unit_num, *) "CaseType:", trim(self % case_type)
    write(unit_num, *) "WindSpeed:", self % WindSpeed
    write(unit_num, *) "WindDir:", self % WindDir
    write(unit_num, *) "Depth:", self % Depth
    write(unit_num, *) "AxisPos:", self % AxisPos
    write(unit_num, *) "BariPos:", self % BariPos
    write(unit_num, *) "In_farm:", self % in_farm
    write(unit_num, *) "Nm:", self % Nmembers
    write(unit_num, *) "Nn:", self % Nnodes
    
    ! Delegate parameter save to the acoustic solver
    if (ALLOCATED(self % acoustic_solver)) then
        write(unit_num, *) "Method:", trim(self % acoustic_solver % get_name())
        CALL self % acoustic_solver % save_parameters(unit_num)
    end if
END SUBROUTINE save_parameters


SUBROUTINE check_acoustic_solver(self)
    CLASS(WindTurbine_t), INTENT(IN) :: self

    LOGICAL :: is_allocated


    if (.not. ALLOCATED(self % acoustic_solver)) then
        error stop "An acoustic solver has to be assigned before pressure computations."
    end if

END SUBROUTINE check_acoustic_solver


SUBROUTINE check_observers_distances(self, observers, min_distance)
    USE Kinds, ONLY: WP, I32
    CLASS(WindTurbine_t), INTENT(IN) :: self
    REAL(WP)            , INTENT(IN) :: observers(:,:)          ! (Nobs, 3)
    REAL(WP)            , INTENT(IN), OPTIONAL :: min_distance  ! [m]

    ! Local variables
    REAL(WP) :: min_dist
    INTEGER(I32) :: i, j, n_nodes, n_obs
    REAL(WP) :: dx, dy, dz, dist

    ! Default min value
    min_dist = 1.0_WP; if (PRESENT(min_distance)) min_dist = min_distance
    
    n_nodes = size(self % x, 1)
    n_obs   = size(observers, 1)
    
    do i = 1, n_obs
        do j = 1, n_nodes
            dx = self % x(j,1) - observers(i,1)
            dy = self % x(j,2) - observers(i,2)
            dz = self % x(j,3) - observers(i,3)
            dist = sqrt(dx*dx + dy*dy + dz*dz)

            if (abs(dist - min_dist)<1e-6_WP) then
                write(*, '(A,I0,A,3(F12.4,1X),A,F6.2,A)') &
                "Observer ", i, " at (", observers(i,:), ") is too close to a turbine node (min distance: ", min_dist, " m)."
                error stop
            end if
        end do
    end do

END SUBROUTINE check_observers_distances


! ---------- RUN FUNCTIONS ---------- !
SUBROUTINE run_spectrums(self, observers, z_obs, block_size)
    CLASS(WindTurbine_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: observers(:,:)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: z_obs
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)                 :: z_obs_
    INTEGER(I32)             :: block_size_

    CALL self % check_acoustic_solver()
    print*, ''; print*, 'Computing spectrums at observer points...'

    ! Defaults
    z_obs_ = - self % Depth/2.0_WP; if (PRESENT(z_obs)) z_obs_ = z_obs
    block_size_ = 4_I32; if (PRESENT(block_size)) block_size_ = block_size
    if (PRESENT(observers)) then
        observers_ = observers
    else
        ALLOCATE(observers_(3,3))
        observers_(:,3) = z_obs_
        observers_(:,2) = 0.0_WP
        observers_(:,1) = [10.0_WP, 250.0_WP, 500.0_WP]
    end if

    ! Compute pressure
    CALL self % check_observers_distances(observers_)
    CALL self % acoustic_solver % compute_pressure(observers_, block_size_, p)

    if (self % debug) print*, "Spectrum Pressure norm: ", norm2(abs(p))

END SUBROUTINE run_spectrums


SUBROUTINE run_polar(self, r, z, n_theta, center, block_size)
    CLASS(WindTurbine_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: r
    REAL(WP)    , INTENT(IN), OPTIONAL  :: z
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_theta
    REAL(WP)    , INTENT(IN), OPTIONAL  :: center(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)                 :: r_, z_, center_(2), theta_deg, theta_rad, d_theta
    INTEGER(I32)             :: n_theta_, block_size_, i

    CALL self % check_acoustic_solver()

    ! Defaults
    r_ = 500.0_WP; if (PRESENT(r)) r_ = r
    z_ = - self % Depth/2.0_WP; if (PRESENT(z)) z_ = z
    n_theta_ = 72_I32; if (PRESENT(n_theta)) n_theta_ = n_theta
    block_size_ = 16_I32; if (PRESENT(block_size)) block_size_ = block_size
        
    center_ = [self % BariPos(1), self % BariPos(2)]
    if (PRESENT(center)) center_ = center

    ! Build observers array
    ALLOCATE(observers_(n_theta_, 3))
    d_theta = 360.0_WP / REAL(n_theta_, WP)
        
    do i = 1, n_theta_
        theta_deg = REAL(i - 1, WP) * d_theta
        theta_rad = theta_deg * PI / 180.0_WP
        observers_(i, 1) = r_ * COS(theta_rad) + center_(1)
        observers_(i, 2) = r_ * SIN(theta_rad) + center_(2)
        observers_(i, 3) = z_
    end do

    CALL self % check_observers_distances(observers_)

    print*, ''
    print*, 'Computing polar pressure (r=', r_, 'm, center=(', center_(1), ','&
            , center_(2), ') m, z=', z_, 'm), observers=', n_theta_, '...'
        
    CALL self % acoustic_solver % compute_pressure(observers_, block_size_, p)
    if (self % debug) print*, "Polar Pressure Norm: ", norm2(abs(p))

END SUBROUTINE run_polar


SUBROUTINE run_cylinder(self, r, n_theta, nz, center, block_size)
    CLASS(WindTurbine_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: r
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_theta
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nz
    REAL(WP)    , INTENT(IN), OPTIONAL  :: center(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

        ! Local variables
     COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)                 :: r_, center_(2), theta_deg, theta_rad, d_theta, z_val, dz
    INTEGER(I32)             :: n_theta_, nz_, block_size_, i, j, idx, Nobs

    CALL self % check_acoustic_solver()

    ! Defaults
    r_ = 500.0_WP; if (PRESENT(r)) r_ = r
    n_theta_ = 72_I32; if (PRESENT(n_theta)) n_theta_ = n_theta
    nz_ = 20_I32; if (PRESENT(nz)) nz_ = nz
    block_size_ = 256_I32; if (PRESENT(block_size)) block_size_ = block_size
        
    center_ = [self % BariPos(1), self % BariPos(2)]
    if (PRESENT(center)) center_ = center

    Nobs = n_theta_ * nz_
    ALLOCATE(observers_(Nobs, 3))
        
    d_theta = 360.0_WP / REAL(n_theta_, WP)
    dz = (0.0_WP - (-self % Depth)) / MAX(REAL(nz_ - 1, WP), 1.0_WP)

    ! Build observers array (meshgrid equivalent)
    idx = 1
    do i = 1, n_theta_
        theta_deg = REAL(i - 1, WP) * d_theta
        theta_rad = theta_deg * PI / 180.0_WP
        do j = 1, nz_
            z_val = -self % Depth + REAL(j - 1, WP) * dz
            observers_(idx, 1) = r_ * COS(theta_rad) + center_(1)
            observers_(idx, 2) = r_ * SIN(theta_rad) + center_(2)
            observers_(idx, 3) = z_val
            idx = idx + 1
        end do
    end do

    CALL self % check_observers_distances(observers_)

    print*, ''
    print*, 'Computing cylindrical pressure (r=', r_, 'm, center=(', center_(1), ','&
            , center_(2), ') m, grid=', n_theta_, 'x', nz_, '=', Nobs, 'obs)...'

    CALL self % acoustic_solver % compute_pressure(observers_, block_size_, p)
    if (self % debug) print*, "Cylinder Pressure Norm: ", norm2(abs(p))

END SUBROUTINE run_cylinder


SUBROUTINE run_decay(self, distance, n_points, z, logspace, block_size)
    CLASS(WindTurbine_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: distance(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_points
    REAL(WP)    , INTENT(IN), OPTIONAL  :: z
    LOGICAL     , INTENT(IN), OPTIONAL  :: logspace
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)                 :: distance_(2), z_, wind_unit(2), d_val, WindDir_rad
    REAL(WP)                 :: d_log1, d_log2, d_step, linear_step
    INTEGER(I32)             :: n_points_, block_size_, i
    LOGICAL                  :: logspace_

    CALL self % check_acoustic_solver()

    ! Defaults
    z_ = -self % Depth/2.0_WP; if (PRESENT(z)) z_ = z
    distance_ = [10.0_WP, 500.0_WP]; if (PRESENT(distance)) distance_ = distance
    n_points_ = 200_I32; if (PRESENT(n_points)) n_points_ = n_points
    logspace_ = .TRUE.; if (PRESENT(logspace)) logspace_ = logspace
    block_size_ = 64_I32; if (PRESENT(block_size)) block_size_ = block_size

    ! Wind direction
    WindDir_rad = self % WindDir * PI / 180.0_WP
    wind_unit = [COS(WindDir_rad), SIN(WindDir_rad)]

    ALLOCATE(observers_(n_points_, 3))

    ! Distances along the line
    if (logspace_) then
        d_log1 = LOG10(distance_(1))
        d_log2 = LOG10(distance_(2))
        d_step = (d_log2 - d_log1) / REAL(n_points_ - 1, WP)
        do i = 1, n_points_
            d_val = 10.0_WP ** (d_log1 + REAL(i - 1, WP) * d_step)
            observers_(i, 1:2) = self % AxisPos(1:2) + d_val * wind_unit
            observers_(i, 3) = z_
        end do
    else
        linear_step = (distance_(2) - distance_(1)) / REAL(n_points_ - 1, WP)
        do i = 1, n_points_
            d_val = distance_(1) + REAL(i - 1, WP) * linear_step
            observers_(i, 1:2) = self % AxisPos(1:2) + d_val * wind_unit
            observers_(i, 3) = z_
        end do
    end if

    CALL self % check_observers_distances(observers_)

    print*, ''
    print*, 'Computing distance decay (points=', n_points_, ', direction=', self % WindDir, 'deg, z=', z_, 'm)...'

    CALL self % acoustic_solver % compute_pressure(observers_, block_size_, p)
    if (self % debug) print*, "Decay Pressure Norm: ", norm2(abs(p))

END SUBROUTINE run_decay


SUBROUTINE run_line(self, p1, p2, n_points, logspace, block_size)
    CLASS(WindTurbine_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: p1(3)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: p2(3)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_points
    LOGICAL     , INTENT(IN), OPTIONAL  :: logspace
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)                 :: p1_(3), p2_(3), dir_vec(3), dist_total, d_val
    REAL(WP)                 :: d_log1, d_log2, d_step, linear_step
    INTEGER(I32)             :: n_points_, block_size_, i
    LOGICAL                  :: logspace_

    CALL self % check_acoustic_solver()

    ! Defaults
    p1_ = [100.0_WP, 0.0_WP, 0.0_WP]; if (PRESENT(p1)) p1_ = p1
    p2_ = [100.0_WP, 0.0_WP, -self % Depth]; if (PRESENT(p2)) p2_ = p2
    logspace_ = .FALSE.; if (PRESENT(logspace)) logspace_ = logspace

    dist_total = NORM2(p2_ - p1_)
    if (dist_total == 0.0_WP) then
        print*, "ERROR: Windturbine.run_line(): points p1 and p2 are coincident."
        STOP
    end if

    if (PRESENT(n_points)) then
        n_points_ = n_points
    else
        if (dist_total > 200.0_WP) then
            n_points_ = INT(dist_total / 5.0_WP, I32)
        else
            n_points_ = INT(dist_total, I32)
        end if
    end if

    if (PRESENT(block_size)) then
        block_size_ = block_size
    else
        block_size_ = INT(n_points_ / 10_I32, I32)
        if (block_size_ < 1_I32) block_size_ = 1_I32
    end if

    ALLOCATE(observers_(n_points_, 3))
    dir_vec = (p2_ - p1_) / dist_total

    ! Distances along the line
    if (logspace_) then
        d_log1 = LOG10(dist_total * 1e-6_WP)
        d_log2 = LOG10(dist_total)
        d_step = (d_log2 - d_log1) / REAL(n_points_ - 1, WP)
        do i = 1, n_points_
            d_val = 10.0_WP ** (d_log1 + REAL(i - 1, WP) * d_step)
            observers_(i, :) = p1_ + (d_val / dist_total) * dir_vec
        end do
    else
        linear_step = dist_total / REAL(n_points_ - 1, WP)
        do i = 1, n_points_
            d_val = REAL(i - 1, WP) * linear_step
            observers_(i, :) = p1_ + (d_val / dist_total) * dir_vec
        end do
    end if

    CALL self % check_observers_distances(observers_)

    print*, ''
    print*, 'Computing line (points=', n_points_, '), p1=', p1_, 'm, p2=', p2_, 'm)...'

    CALL self % acoustic_solver % compute_pressure(observers_, block_size_, p)
    if (self % debug) print*, "Line Pressure Norm: ", norm2(abs(p))

END SUBROUTINE run_line


SUBROUTINE run_sliceXY(self, z, nx, ny, xlim, ylim, center, block_size)
    CLASS(WindTurbine_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: z
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nx
    INTEGER(I32), INTENT(IN), OPTIONAL  :: ny
    REAL(WP)    , INTENT(IN), OPTIONAL  :: xlim(2)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: ylim(2)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: center(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)                 :: z_, xlim_(2), ylim_(2), center_(2)
    REAL(WP)                 :: mid_x, mid_y, dx, dy, x_val, y_val
    INTEGER(I32)             :: nx_, ny_, block_size_, i, j, idx
    REAL(WP), PARAMETER      :: TOL = 1e-5_WP

    CALL self % check_acoustic_solver()

    ! Defaults
    z_ = -self % Depth/2.0_WP; if (PRESENT(z)) z_ = z
    if (z_ > 0.0_WP) then
        print*, "WindTurbine.run_sliceXY(): z must be <= 0. Switching to -z"
        z_ = -z_
    end if

    nx_ = 26_I32; if (PRESENT(nx)) nx_ = nx
    ny_ = 26_I32; if (PRESENT(ny)) ny_ = ny
    xlim_ = [-500.0_WP, 500.0_WP]; if (PRESENT(xlim)) xlim_ = xlim
    ylim_ = [-500.0_WP, 500.0_WP]; if (PRESENT(ylim)) ylim_ = ylim
    block_size_ = 128_I32; if (PRESENT(block_size)) block_size_ = block_size

    mid_x = 0.5_WP * (xlim_(1) + xlim_(2))
    mid_y = 0.5_WP * (ylim_(1) + ylim_(2))

    if (PRESENT(center)) then
        if (ABS(center(1) - mid_x) > TOL .OR. ABS(center(2) - mid_y) > TOL) then
            print*, "WindTurbine.run_sliceXY(): provided center does not match midpoint. Switching to proper center."
            center_ = [mid_x, mid_y]
        else
            center_ = center
        end if
    else
        center_ = [mid_x, mid_y]
    end if

    ALLOCATE(observers_(nx_ * ny_, 3))
    dx = (xlim_(2) - xlim_(1)) / MAX(REAL(nx_ - 1, WP), 1.0_WP)
    dy = (ylim_(2) - ylim_(1)) / MAX(REAL(ny_ - 1, WP), 1.0_WP)

    ! Build observers (meshgrid equivalent)
    idx = 1
    do i = 1, ny_
        y_val = ylim_(1) + REAL(i - 1, WP) * dy
        do j = 1, nx_
            x_val = xlim_(1) + REAL(j - 1, WP) * dx
            observers_(idx, 1) = x_val
            observers_(idx, 2) = y_val
            observers_(idx, 3) = z_
            idx = idx + 1
        end do
    end do

    CALL self % check_observers_distances(observers_, 0.1_WP)

    print*, ''
    print*, 'Computing XY slice (nx=', nx_, ', ny=', ny_, ', z=', z_, 'm, center=(', center_(1), ',', center_(2), ') m)...'

    CALL self % acoustic_solver % compute_pressure(observers_, block_size_, p)
    if (self % debug) print*, "SliceXY Pressure Norm: ", norm2(abs(p))

END SUBROUTINE run_sliceXY


SUBROUTINE run_sliceXZ(self, y, nx, nz, xlim, zlim, block_size)
    CLASS(WindTurbine_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: y
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nx
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nz
    REAL(WP)    , INTENT(IN), OPTIONAL  :: xlim(2)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: zlim(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)                 :: y_, xlim_(2), zlim_(2), dx, dz, x_val, z_val
    INTEGER(I32)             :: nx_, nz_, block_size_, i, j, idx

    CALL self % check_acoustic_solver()

    ! Defaults
    y_ = self % AxisPos(2); if (PRESENT(y)) y_ = y
    nx_ = 26_I32; if (PRESENT(nx)) nx_ = nx
    nz_ = 26_I32; if (PRESENT(nz)) nz_ = nz
    xlim_ = [-500.0_WP, 500.0_WP]; if (PRESENT(xlim)) xlim_ = xlim
    zlim_ = [0.0_WP, -self % Depth]; if (PRESENT(zlim)) zlim_ = zlim
    block_size_ = 128_I32; if (PRESENT(block_size)) block_size_ = block_size

    ALLOCATE(observers_(nx_ * nz_, 3))
    dx = (xlim_(2) - xlim_(1)) / MAX(REAL(nx_ - 1, WP), 1.0_WP)
    dz = (zlim_(2) - zlim_(1)) / MAX(REAL(nz_ - 1, WP), 1.0_WP)

    ! Build observers
    idx = 1
    do i = 1, nz_
        z_val = zlim_(1) + REAL(i - 1, WP) * dz
        do j = 1, nx_
            x_val = xlim_(1) + REAL(j - 1, WP) * dx
            observers_(idx, 1) = x_val
            observers_(idx, 2) = y_
            observers_(idx, 3) = z_val
            idx = idx + 1
        end do
    end do

    CALL self % check_observers_distances(observers_, 0.1_WP)

    print*, ''
    print*, 'Computing XZ slice (nx=', nx_, ', nz=', nz_, ', y=', y_, 'm, z from', zlim_(1), 'to', zlim_(2), 'm)...'

    CALL self % acoustic_solver % compute_pressure(observers_, block_size_, p)
    if (self % debug) print*, "SliceXZ Pressure Norm: ", norm2(abs(p))

END SUBROUTINE run_sliceXZ


! SUBROUTINE run_sphere(self, r, n_theta, nz, center, block_size)
!     USE MethodImages, ONLY: MethodImages_t

!     CLASS(WindTurbine_t), INTENT(INOUT) :: self
!     REAL(WP)    , INTENT(IN), OPTIONAL  :: r
!     INTEGER(I32), INTENT(IN), OPTIONAL  :: n_theta
!     INTEGER(I32), INTENT(IN), OPTIONAL  :: nz
!     REAL(WP)    , INTENT(IN), OPTIONAL  :: center(3)
!     INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

!     ! Local variables
!     COMPLEX(WP), ALLOCATABLE :: p(:,:)
!     REAL(WP)   , ALLOCATABLE :: observers_(:,:)
!     REAL(WP)                 :: r_, center_(3), max_dist, current_dist
!     REAL(WP)                 :: dz, dtheta, zc, r_xy, theta_c
!     INTEGER(I32)             :: n_theta_, nz_, block_size_, i, j, idx, Nobs

!     CALL self % check_acoustic_solver()

!     SELECT TYPE (solver => self % acoustic_solver)
!     TYPE IS (MethodImages_t)

!         ! Defaults
!         r_ = 30.0_WP
!         if (PRESENT(r)) r_ = r
!         n_theta_ = 72_I32
!         if (PRESENT(n_theta)) n_theta_ = n_theta
!         nz_ = 20_I32
!         if (PRESENT(nz)) nz_ = nz
!         block_size_ = 128_I32
!         if (PRESENT(block_size)) block_size_ = block_size

!         if (PRESENT(center)) then
!             center_ = center
!         else
!             center_ = [self % BariPos(1), self % BariPos(2), -self % Depth / 2.0_WP]
!         end if

!         ! Basic enclosure check
!         max_dist = 0.0_WP
!         do i = 1, SIZE(self % x, 1)
!             current_dist = NORM2(self % x(i,:) - center_)
!             if (current_dist > max_dist) max_dist = current_dist
!         end do
!         max_dist = max_dist * 1.1_WP

!         if (r_ < max_dist) then
!             print*, "WindTurbine.run_sphere(): requested radius r=", r_, "m is smaller than max dist."
!             print*, "Radius increased to ", max_dist * 1.05_WP, "m"
!             r_ = max_dist * 1.05_WP
!         end if

!         Nobs = n_theta_ * nz_
!         ALLOCATE(observers_(Nobs, 3))

!         dz = (2.0_WP * r_) / REAL(nz_, WP)
!         dtheta = (2.0_WP * PI) / REAL(n_theta_, WP)

!         idx = 1
!         do i = 1, nz_
!             zc = (center_(3) - r_) + REAL(i - 1, WP) * dz + (dz / 2.0_WP)
!             r_xy = SQRT(MAX(0.0_WP, r_**2 - (zc - center_(3))**2))

!             do j = 1, n_theta_
!                 theta_c = REAL(j - 1, WP) * dtheta + (dtheta / 2.0_WP)
!                 observers_(idx, 1) = center_(1) + r_xy * COS(theta_c)
!                 observers_(idx, 2) = center_(2) + r_xy * SIN(theta_c)
!                 observers_(idx, 3) = zc
!                 idx = idx + 1
!             end do
!         end do

!         CALL self % check_observers_distances(observers_, 0.1_WP)

!         print*, ''
!         print*, 'Computing spherical pressure field ...'
!         print*, '  Centre: (', center_(1), ',', center_(2), ',', center_(3), ') m'
!         print*, '  Radius: ', r_, 'm | grid:', n_theta_, 'azim x', nz_, 'z (', Nobs, 'observers)'

!         CALL solver % set_N_images(0_I32)
!         CALL solver % compute_pressure(observers_, block_size_, p)
!         CALL solver % restore_default_images()

!         if (self % debug) print*, "Sphere Pressure Norm: ", norm2(abs(p))

!     CLASS DEFAULT
!         ! This geometry is only implemented for MethodImages_t.
!         RETURN
!     END SELECT

! END SUBROUTINE run_sphere


SUBROUTINE run_all(self, &
    spectrums_observers, spectrums_z_obs, spectrums_block_size, &
    polar_r, polar_z, polar_n_theta, polar_center, polar_block_size, &
    cylinder_r, cylinder_n_theta, cylinder_nz, cylinder_center, cylinder_block_size, &
    decay_distance, decay_n_points, decay_z, decay_logspace, decay_block_size, &
    line_p1, line_p2, line_n_points, line_logspace, line_block_size, &
    sliceXY_z, sliceXY_nx, sliceXY_ny, sliceXY_xlim, sliceXY_ylim, sliceXY_center, sliceXY_block_size, &
    sliceXZ_y, sliceXZ_nx, sliceXZ_nz, sliceXZ_xlim, sliceXZ_zlim, sliceXZ_block_size, &
    sphere_r, sphere_n_theta, sphere_nz, sphere_center, sphere_block_size)
        
    CLASS(WindTurbine_t), INTENT(INOUT) :: self
    ! --- spectrums ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: spectrums_observers(:,:)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: spectrums_z_obs
    INTEGER(I32), INTENT(IN), OPTIONAL  :: spectrums_block_size
    ! --- polar ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: polar_r, polar_z, polar_center(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: polar_n_theta, polar_block_size
    ! --- cylinder ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: cylinder_r, cylinder_center(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: cylinder_n_theta, cylinder_nz, cylinder_block_size
    ! --- decay ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: decay_distance(2), decay_z
    INTEGER(I32), INTENT(IN), OPTIONAL  :: decay_n_points, decay_block_size
    LOGICAL     , INTENT(IN), OPTIONAL  :: decay_logspace
    ! --- line ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: line_p1(3), line_p2(3)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: line_n_points, line_block_size
    LOGICAL     , INTENT(IN), OPTIONAL  :: line_logspace
    ! --- sliceXY ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: sliceXY_z, sliceXY_xlim(2), sliceXY_ylim(2), sliceXY_center(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: sliceXY_nx, sliceXY_ny, sliceXY_block_size
    ! --- sliceXZ ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: sliceXZ_y, sliceXZ_xlim(2), sliceXZ_zlim(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: sliceXZ_nx, sliceXZ_nz, sliceXZ_block_size
    ! --- sphere ---
    REAL(WP)    , INTENT(IN), OPTIONAL  :: sphere_r, sphere_center(3)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: sphere_n_theta, sphere_nz, sphere_block_size

    CALL self % check_acoustic_solver()

    print*, ''
    print*, '=================================================='
    print*, ' RUNNING ALL ACOUSTIC POST-PROCESSING'
    print*, ' Total steps: 8'
    print*, '=================================================='

    print*, ''
    print*, '[Step 1/8] Spectrums at observer points ...'
    CALL self % run_spectrums(spectrums_observers, spectrums_z_obs, spectrums_block_size)

    print*, ''
    print*, '[Step 2/8] Polar contour ...'
    CALL self % run_polar(polar_r, polar_z, polar_n_theta, polar_center, polar_block_size)

    print*, ''
    print*, '[Step 3/8] Distance decay along wind ...'
    CALL self % run_decay(decay_distance, decay_n_points, decay_z, decay_logspace, decay_block_size)

    print*, ''
    print*, '[Step 4/8] Line between two points ...'
    CALL self % run_line(line_p1, line_p2, line_n_points, line_logspace, line_block_size)

    print*, ''
    print*, '[Step 5/8] Spherical radiation (free-field) ...'
    ! CALL self % run_sphere(sphere_r, sphere_n_theta, sphere_nz, sphere_center, sphere_block_size)

    print*, ''
    print*, '[Step 6/8] Horizontal slice XY ...'
    CALL self % run_sliceXY(sliceXY_z, sliceXY_nx, sliceXY_ny, sliceXY_xlim, sliceXY_ylim, sliceXY_center, sliceXY_block_size)

    print*, ''
    print*, '[Step 7/8] Vertical slice XZ ...'
    CALL self % run_sliceXZ(sliceXZ_y, sliceXZ_nx, sliceXZ_nz, sliceXZ_xlim, sliceXZ_zlim, sliceXZ_block_size)

    print*, ''
    print*, '[Step 8/8] Cylinder surface ...'
    CALL self % run_cylinder(cylinder_r, cylinder_n_theta, cylinder_nz, cylinder_center, cylinder_block_size)

    print*, ''
    print*, '=================================================='
    print*, ' ALL ACOUSTIC COMPUTATIONS FINISHED'
    print*, '=================================================='

END SUBROUTINE run_all


END MODULE WindTurbine