MODULE AcousticSolver

USE Kinds, ONLY: WP, I32, PI, SPEED_OF_SOUND, RHO_WATER, PREF => P_REF, IN_CLUSTER
USE IOUtils, ONLY: save_to_hdf5
USE WindTurbine, ONLY: WindTurbine_t

IMPLICIT NONE

PRIVATE
PUBLIC :: AcousticSolver_t, save_turbine_params


! ------------------
! Abstract Base Type
! ------------------
TYPE, ABSTRACT :: AcousticSolver_t
    
    REAL(WP)          :: c       = SPEED_OF_SOUND       ! [m/s] Speed of sound in fluid
    REAL(WP)          :: rho     = RHO_WATER            ! [kg/m^3] Density of the fluid
    REAL(WP)          :: eps     = 1e-12_WP             ! [-] Small parameter for regularization
    REAL(WP)          :: p_ref   = PREF                 ! [Pa] Reference pressure
    LOGICAL           :: cluster = IN_CLUSTER           ! [-] Flag for cluster execution
    LOGICAL           :: debug   = .false.              ! [-] Debug flag
    CHARACTER(len=50) :: solver_name = 'None'           ! [-] Solver name
    CHARACTER(len=512):: save_path = 'None'             ! [-] Path to save acoustics

    ! Wind Turbine System PlaceHolder
    CLASS(WindTurbine_t), ALLOCATABLE :: turbines(:)

    CONTAINS

    ! --- Public type-bound procedures --- !
    PROCEDURE :: save_turbine_params            ! Saves turbine parameters to save_path
    PROCEDURE :: save_acoustics
    PROCEDURE :: check_observers_distances      ! Helper for run subroutines
    PROCEDURE :: run_spectrums                  ! Subroutines to compute acoustic pressure
    PROCEDURE :: run_polar
    PROCEDURE :: run_cylinder
    PROCEDURE :: run_decay
    PROCEDURE :: run_line
    PROCEDURE :: run_sliceXY
    PROCEDURE :: run_sliceXZ
    PROCEDURE :: run_all
    PROCEDURE :: run_sphere => dummy_run_sphere

    ! --- Deferred(abstract) procedures --- !
    PROCEDURE(get_name_i)        , DEFERRED :: get_name
    PROCEDURE(compute_pressure_i), DEFERRED :: compute_pressure
    PROCEDURE(save_parameters_i) , DEFERRED :: save_parameters

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

! Save parameters
ABSTRACT INTERFACE
    SUBROUTINE save_parameters_i(self)
        IMPORT:: AcousticSolver_t
        CLASS(AcousticSolver_t), INTENT(IN) :: self
    END SUBROUTINE
END INTERFACE


CONTAINS

! ---------- HELPERS ---------- !
SUBROUTINE save_turbine_params(self)
    CLASS(AcousticSolver_t), INTENT(IN) :: self

    ! Local variables
    CHARACTER(len=512) :: save_path
    INTEGER(I32)       :: Nturb, i
    CHARACTER(len=2)   :: ichar

    save_path = trim(self % save_path)
    Nturb     = size(self % turbines)

    ! Turbines
    if (Nturb > 1) then
        do i = 1, Nturb
            WRITE(ichar, '(I2.2)') i
            CALL save_to_hdf5(save_path, "WindSpeed_"       // ichar, self % turbines(i) % WindSpeed)
            CALL save_to_hdf5(save_path, "WindDir_"         // ichar, self % turbines(i) % WindDir)
            CALL save_to_hdf5(save_path, "Depth_"           // ichar, self % turbines(i) % Depth)
            CALL save_to_hdf5(save_path, "Structure_nodes_" // ichar, self % turbines(i) % x_all)
            CALL save_to_hdf5(save_path, "Case_type_"       // ichar, self % turbines(i) % case_type)
            CALL save_to_hdf5(save_path, "In_farm_"         // ichar, self % turbines(i) % in_farm)
            CALL save_to_hdf5(save_path, "AxisPos_"         // ichar, self % turbines(i) % AxisPos)
            CALL save_to_hdf5(save_path, "BariPos_"         // ichar, self % turbines(i) % BariPos)
            CALL self % turbines(i) % save_specific_params(save_path, ichar)
        end do
    elseif (Nturb == 1) then
        CALL save_to_hdf5(save_path, "WindSpeed"       , self % turbines(1) % WindSpeed)
        CALL save_to_hdf5(save_path, "WindDir"         , self % turbines(1) % WindDir)
        CALL save_to_hdf5(save_path, "Depth"           , self % turbines(1) % Depth)
        CALL save_to_hdf5(save_path, "Structure_nodes" , self % turbines(1) % x_all)
        CALL save_to_hdf5(save_path, "Case_type"       , self % turbines(1) % case_type)
        CALL save_to_hdf5(save_path, "In_farm"         , self % turbines(1) % in_farm)
        CALL save_to_hdf5(save_path, "AxisPos"         , self % turbines(1) % AxisPos)
        CALL save_to_hdf5(save_path, "BariPos"         , self % turbines(1) % BariPos)
        CALL self % turbines(1) % save_specific_params(save_path)
    end if

END SUBROUTINE save_turbine_params


SUBROUTINE save_acoustics(Self, var_name, var)
    CLASS(AcousticSolver_t), INTENT(IN) :: self
    CHARACTER(len=*)       , INTENT(IN) :: var_name
    CLASS(*)               , INTENT(IN) :: var(..)

    ! Local variables
    CHARACTER(len=512) :: save_path

    save_path = trim(self % save_path)
    
    SELECT RANK (v_rank => var)
    RANK(0) 
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP)) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP)); CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (LOGICAL) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        END SELECT
    RANK(1) 
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP)) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP)); CALL save_to_hdf5(save_path, trim(var_name), v)
        END SELECT
    RANK(2) 
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP)) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP)); CALL save_to_hdf5(save_path, trim(var_name), v)
        END SELECT
    RANK(3) 
        SELECT TYPE (v => v_rank)
        TYPE IS (REAL(WP)) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (INTEGER) ; CALL save_to_hdf5(save_path, trim(var_name), v)
        TYPE IS (COMPLEX(WP)); CALL save_to_hdf5(save_path, trim(var_name), v)
        END SELECT
    END SELECT

END SUBROUTINE save_acoustics


SUBROUTINE check_observers_distances(self, observers, min_distance)
    CLASS(AcousticSolver_t), INTENT(IN) :: self
    REAL(WP)               , INTENT(IN) :: observers(:,:)          ! (Nobs, 3)
    REAL(WP)               , INTENT(IN), OPTIONAL :: min_distance  ! [m]

    ! Local variables
    REAL(WP) :: min_dist
    INTEGER(I32) :: i, j, n_nodes, n_obs
    REAL(WP) :: dx, dy, dz, dist

    ! Default min value
    min_dist = 1.0_WP; if (PRESENT(min_distance)) min_dist = min_distance
    
    n_nodes = size(self % turbines(1) % x, 1)
    n_obs   = size(observers, 1)
    
    do i = 1, n_obs
        do j = 1, n_nodes
            dx = self % turbines(1) % x(j,1) - observers(i,1)
            dy = self % turbines(1) % x(j,2) - observers(i,2)
            dz = self % turbines(1) % x(j,3) - observers(i,3)
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
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: observers(:,:)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: z_obs
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)                 :: z_obs_
    INTEGER(I32)             :: block_size_

    print*, ''; print*, 'Computing spectrums at observer points...'

    ! Defaults
    z_obs_ = - self % turbines(1) % Depth/2.0_WP; if (PRESENT(z_obs))      z_obs_      = z_obs
    block_size_ = 4_I32                         ; if (PRESENT(block_size)) block_size_ = block_size
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
    CALL self % compute_pressure(observers_, block_size_, p)

    if (self % debug) write(*, '(A, F0.2)') 'Spectrum Pressure norm: ', norm2(abs(p))

    ! Save data
    CALL self % save_acoustics("P_spectrums", p)
    CALL self % save_acoustics("Obs_spectrums", observers_)
    CALL self % save_acoustics("Freqs", self % turbines(1) % Freqs)

    print*, ''
    print*, '   ---> Spectrum data saved at ', self % save_path

END SUBROUTINE run_spectrums


SUBROUTINE run_polar(self, r, z, n_theta, center, block_size)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: r
    REAL(WP)    , INTENT(IN), OPTIONAL  :: z
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_theta
    REAL(WP)    , INTENT(IN), OPTIONAL  :: center(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:), theta_deg_polar(:)
    REAL(WP)                 :: r_, z_, center_(2), theta_deg, theta_rad, d_theta
    INTEGER(I32)             :: n_theta_, block_size_, i

    ! Defaults
    r_ = 500.0_WP                           ; if (PRESENT(r))          r_          = r
    z_ = - self % turbines(1) % Depth/2.0_WP; if (PRESENT(z))          z_          = z
    n_theta_ = 72_I32                       ; if (PRESENT(n_theta))    n_theta_    = n_theta
    block_size_ = 16_I32                    ; if (PRESENT(block_size)) block_size_ = block_size
        
    center_ = [self % turbines % BariPos(1), self % turbines % BariPos(2)]
    if (PRESENT(center)) center_ = center

    ! Build observers array
    ALLOCATE(observers_(n_theta_, 3), theta_deg_polar(n_theta_))
    d_theta = 360.0_WP / REAL(n_theta_, WP)
        
    do i = 1, n_theta_
        theta_deg = REAL(i - 1, WP) * d_theta
        theta_rad = theta_deg * PI / 180.0_WP
        observers_(i, 1) = r_ * COS(theta_rad) + center_(1)
        observers_(i, 2) = r_ * SIN(theta_rad) + center_(2)
        observers_(i, 3) = z_
        theta_deg_polar(i) = theta_deg
    end do

    CALL self % check_observers_distances(observers_)

    write(*, '(A)') ''
    write(*, '(A, F0.2, A, F0.2, A, F0.2, A, F0.2, A, I0, A)') &
        'Computing polar pressure (r=', r_, ' m, center=(', center_(1), ', ', center_(2), &
        ') m, z=', z_, ' m), observers=', n_theta_, '...'
        
    ! Compute_pressure
    CALL self % compute_pressure(observers_, block_size_, p)
    if (self % debug) write(*, '(A, F0.2)') 'Polar Pressure Norm: ', norm2(abs(p))

    ! Save data
    CALL self % save_acoustics("Freqs", self % turbines(1) % Freqs)
    CALL self % save_acoustics("P_polar", p)
    CALL self % save_acoustics("R_polar", r_)
    CALL self % save_acoustics("Z_polar", z_)
    CALL self % save_acoustics("Theta_deg_polar", theta_deg_polar)
    CALL self % save_acoustics("Obs_polar", observers_)
    CALL self % save_acoustics("Center_polar", center_)

    print*, ''
    print*, '   ---> Polar data saved at ', self % save_path

END SUBROUTINE run_polar


SUBROUTINE run_cylinder(self, r, n_theta, nz, center, block_size)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: r
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_theta
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nz
    REAL(WP)    , INTENT(IN), OPTIONAL  :: center(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)                 :: r_, center_(2), theta_rad, d_theta, z_val, dz
    INTEGER(I32)             :: n_theta_, nz_, block_size_, i, j, idx, Nobs

    ! Additional for saving
    REAL(WP), ALLOCATABLE    :: theta_deg_arr(:)      ! Azimuth angles in degrees
    REAL(WP), ALLOCATABLE    :: z_cylinder(:)         ! Vertical coordinates
    REAL(WP), ALLOCATABLE    :: obs_reshaped(:,:,:)   ! Reshaped observers (n_theta, nz, 3)
    COMPLEX(WP), ALLOCATABLE :: p_reshaped(:,:,:)     ! Reshaped pressure (Nfreqs, n_theta, nz)
    REAL(WP)                 :: dA                    ! Differential area element

    ! Defaults
    r_ = 500.0_WP        ; if (PRESENT(r))          r_          = r
    n_theta_ = 72_I32    ; if (PRESENT(n_theta))    n_theta_    = n_theta
    nz_ = 20_I32         ; if (PRESENT(nz))         nz_         = nz
    block_size_ = 256_I32; if (PRESENT(block_size)) block_size_ = block_size

    center_ = [self % turbines(1) % BariPos(1), self % turbines(1) % BariPos(2)]
    if (PRESENT(center)) center_ = center

    Nobs = n_theta_ * nz_
    ALLOCATE(observers_(Nobs, 3))

    d_theta = 360.0_WP / REAL(n_theta_, WP)
    dz = (0.0_WP - (-self % turbines(1) % Depth)) / MAX(REAL(nz_ - 1, WP), 1.0_WP)

    ! Allocate arrays for storing angles and vertical coordinates
    ALLOCATE(theta_deg_arr(n_theta_))
    ALLOCATE(z_cylinder(nz_))

    ! Fill vertical coordinate array
    do j = 1, nz_
        z_cylinder(j) = -self % turbines(1) % Depth + REAL(j - 1, WP) * dz
    end do

    ! Build observers array and angle array (meshgrid equivalent)
    idx = 1
    do i = 1, n_theta_
        theta_deg_arr(i) = REAL(i - 1, WP) * d_theta
        theta_rad = theta_deg_arr(i) * PI / 180.0_WP
        do j = 1, nz_
            z_val = z_cylinder(j)
            observers_(idx, 1) = r_ * COS(theta_rad) + center_(1)
            observers_(idx, 2) = r_ * SIN(theta_rad) + center_(2)
            observers_(idx, 3) = z_val
            idx = idx + 1
        end do
    end do

    CALL self % check_observers_distances(observers_)

    write(*, '(A)') ''
    write(*, '(A, F0.2, A, F0.2, A, F0.2, A, I0, A, I0, A, I0, A)') &
        'Computing cylindrical pressure (r=', r_, ' m, center=(', center_(1), ', ', center_(2), &
        ') m, grid=', n_theta_, ' x ', nz_, ' = ', Nobs, ' obs)...'

    CALL self % compute_pressure(observers_, block_size_, p)
    if (self % debug) write(*, '(A, F0.2)') 'Cylinder Pressure Norm: ', norm2(abs(p))

    ! Reshape pressure and observers to 3D for storage
    ALLOCATE(p_reshaped(SIZE(p,1), n_theta_, nz_))
    ALLOCATE(obs_reshaped(n_theta_, nz_, 3))

    idx = 1
    do i = 1, n_theta_
        do j = 1, nz_
            obs_reshaped(i,j,:) = observers_(idx,:)
            p_reshaped(:,i,j)   = p(:,idx)
            idx = idx + 1
        end do
    end do

    ! Compute differential area element (as in Python)
    dA = r_ * dz * (2.0_WP * PI / REAL(n_theta_, WP))

    ! Save results
    CALL self % save_acoustics("Freqs", self % turbines(1) % Freqs)
    CALL self % save_acoustics("P_cylinder", p_reshaped)
    CALL self % save_acoustics("R_cylinder", r_)
    CALL self % save_acoustics("Z_cylinder", z_cylinder)
    CALL self % save_acoustics("Theta_deg_cylinder", theta_deg_arr)
    CALL self % save_acoustics("Obs_cylinder", obs_reshaped)
    CALL self % save_acoustics("Center_cylinder", center_)
    CALL self % save_acoustics("dA_cylinder", dA)

    print*, ''
    print*, '   ---> Cylinder data saved at ', self % save_path

END SUBROUTINE run_cylinder


SUBROUTINE run_decay(self, distance, n_points, z, logspace, block_size)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: distance(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_points
    REAL(WP)    , INTENT(IN), OPTIONAL  :: z
    LOGICAL     , INTENT(IN), OPTIONAL  :: logspace
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)   , ALLOCATABLE :: distances(:)        ! Actual distances along the line
    REAL(WP)                 :: distance_(2), z_, wind_unit(2), d_val, WindDir_rad
    REAL(WP)                 :: d_log1, d_log2, d_step, linear_step
    INTEGER(I32)             :: n_points_, block_size_, i
    LOGICAL                  :: logspace_
    INTEGER                  :: logspace_int        ! Integer representation for saving

    ! Defaults
    z_ = -self % turbines(1) % Depth/2.0_WP; if (PRESENT(z))          z_          = z
    distance_ = [10.0_WP, 500.0_WP]        ; if (PRESENT(distance))   distance_   = distance
    n_points_ = 200_I32                    ; if (PRESENT(n_points))   n_points_   = n_points
    logspace_ = .true.                     ; if (PRESENT(logspace))   logspace_   = logspace
    block_size_ = 64_I32                   ; if (PRESENT(block_size)) block_size_ = block_size

    ! Wind direction
    WindDir_rad = self % turbines(1) % WindDir * PI / 180.0_WP
    wind_unit = [COS(WindDir_rad), SIN(WindDir_rad)]

    ALLOCATE(observers_(n_points_, 3))
    ALLOCATE(distances(n_points_))

    ! Distances along the line
    if (logspace_) then
        d_log1 = LOG10(distance_(1))
        d_log2 = LOG10(distance_(2))
        d_step = (d_log2 - d_log1) / REAL(n_points_ - 1, WP)
        do i = 1, n_points_
            d_val = 10.0_WP ** (d_log1 + REAL(i - 1, WP) * d_step)
            distances(i) = d_val
            observers_(i, 1:2) = self % turbines(1) % AxisPos(1:2) + d_val * wind_unit
            observers_(i, 3) = z_
        end do
    else
        linear_step = (distance_(2) - distance_(1)) / REAL(n_points_ - 1, WP)
        do i = 1, n_points_
            d_val = distance_(1) + REAL(i - 1, WP) * linear_step
            distances(i) = d_val
            observers_(i, 1:2) = self % turbines(1) % AxisPos(1:2) + d_val * wind_unit
            observers_(i, 3) = z_
        end do
    end if

    CALL self % check_observers_distances(observers_)

    write(*, '(A)') ''
    write(*, '(A, I0, A, F0.2, A, F0.2, A)') &
        'Computing distance decay (points=', n_points_, ', direction=', &
        self % turbines(1) % WindDir, ' deg, z=', z_, ' m)...'

    CALL self % compute_pressure(observers_, block_size_, p)
    if (self % debug) write(*, '(A, F0.2)') 'Decay Pressure Norm: ', norm2(abs(p))

    ! Convert logical to integer for saving
    logspace_int = MERGE(1, 0, logspace_)

    ! Save results
    CALL self % save_acoustics("Freqs", self % turbines(1) % Freqs)
    CALL self % save_acoustics("P_decay", p)
    CALL self % save_acoustics("Z_decay", z_)
    CALL self % save_acoustics("Distance_decay", distance_)
    CALL self % save_acoustics("Distances_decay", distances)
    CALL self % save_acoustics("Obs_decay", observers_)
    CALL self % save_acoustics("Logspace_decay", logspace_int)

    print*, ''
    print*, '   ---> Decay data saved at ', self % save_path

END SUBROUTINE run_decay


SUBROUTINE run_line(self, p1, p2, n_points, logspace, block_size)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: p1(3)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: p2(3)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_points
    LOGICAL     , INTENT(IN), OPTIONAL  :: logspace
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)   , ALLOCATABLE :: distances(:)     ! Distances from p1 along the line
    REAL(WP)                 :: p1_(3), p2_(3), dir_vec(3), dist_total, d_val
    REAL(WP)                 :: d_log1, d_log2, d_step, linear_step
    INTEGER(I32)             :: n_points_, block_size_, i
    LOGICAL                  :: logspace_
    INTEGER                  :: logspace_int     ! Integer representation for saving

    ! Defaults
    p1_ = [100.0_WP, 0.0_WP, 0.0_WP]                     ; if (PRESENT(p1)) p1_ = p1
    p2_ = [100.0_WP, 0.0_WP, -self % turbines(1) % Depth]; if (PRESENT(p2)) p2_ = p2
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
    ALLOCATE(distances(n_points_))
    dir_vec = (p2_ - p1_) / dist_total

    ! Distances along the line
    if (logspace_) then
        d_log1 = LOG10(dist_total * 1e-6_WP)
        d_log2 = LOG10(dist_total)
        d_step = (d_log2 - d_log1) / REAL(n_points_ - 1, WP)
        do i = 1, n_points_
            d_val = 10.0_WP ** (d_log1 + REAL(i - 1, WP) * d_step)
            distances(i) = d_val
            observers_(i, :) = p1_ + (d_val / dist_total) * dir_vec
        end do
    else
        linear_step = dist_total / REAL(n_points_ - 1, WP)
        do i = 1, n_points_
            d_val = REAL(i - 1, WP) * linear_step
            distances(i) = d_val
            observers_(i, :) = p1_ + (d_val / dist_total) * dir_vec
        end do
    end if

    CALL self % check_observers_distances(observers_)

    write(*, '(A)') ''
    write(*, '(A, I0, A, F0.2, A, F0.2, A, F0.2, A, F0.2, A, F0.2, A, F0.2, A)') &
        'Computing line (points=', n_points_, '), p1=(', p1_(1), ', ', p1_(2), ', ', p1_(3), &
        ') m, p2=(', p2_(1), ', ', p2_(2), ', ', p2_(3), ') m)...'

    CALL self % compute_pressure(observers_, block_size_, p)
    if (self % debug) write(*, '(A, F0.2)') 'Line Pressure Norm: ', norm2(abs(p))

    ! Convert logical to integer for saving
    logspace_int = MERGE(1, 0, logspace_)

    ! Save results
    CALL self % save_acoustics("Freqs", self % turbines(1) % Freqs)
    CALL self % save_acoustics("P_line", p)
    CALL self % save_acoustics("P1_line", p1_)
    CALL self % save_acoustics("P2_line", p2_)
    CALL self % save_acoustics("Obs_line", observers_)
    CALL self % save_acoustics("Distance_line", dist_total)
    CALL self % save_acoustics("Distances_line", distances)
    CALL self % save_acoustics("Logspace_line", logspace_int)

    print*, ''
    print*, '   ---> Line data saved at ', self % save_path

END SUBROUTINE run_line


SUBROUTINE run_sliceXY(self, z, nx, ny, xlim, ylim, center, block_size)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: z
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nx
    INTEGER(I32), INTENT(IN), OPTIONAL  :: ny
    REAL(WP)    , INTENT(IN), OPTIONAL  :: xlim(2)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: ylim(2)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: center(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    COMPLEX(WP), ALLOCATABLE :: p_reshaped(:,:,:)   ! (Nfreqs, nx, ny)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)   , ALLOCATABLE :: obs_reshaped(:,:,:) ! (nx, ny, 3)
    REAL(WP)   , ALLOCATABLE :: xs(:), ys(:)        ! Coordinate arrays
    REAL(WP)                 :: z_, xlim_(2), ylim_(2), center_(2)
    REAL(WP)                 :: mid_x, mid_y, dx, dy, x_val, y_val
    INTEGER(I32)             :: nx_, ny_, block_size_, i, j, idx
    REAL(WP), PARAMETER      :: TOL = 1e-5_WP

    ! Defaults
    z_ = -self % turbines(1) % Depth/2.0_WP; if (PRESENT(z)) z_ = z
    if (z_ > 0.0_WP) then
        print*, "WindTurbine.run_sliceXY(): z must be <= 0. Switching to -z"
        z_ = -z_
    end if

    nx_ = 26_I32                 ; if (PRESENT(nx))         nx_         = nx
    ny_ = 26_I32                 ; if (PRESENT(ny))         ny_         = ny
    xlim_ = [-500.0_WP, 500.0_WP]; if (PRESENT(xlim))       xlim_       = xlim
    ylim_ = [-500.0_WP, 500.0_WP]; if (PRESENT(ylim))       ylim_       = ylim
    block_size_ = 128_I32        ; if (PRESENT(block_size)) block_size_ = block_size

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
    ALLOCATE(xs(nx_), ys(ny_))
    dx = (xlim_(2) - xlim_(1)) / MAX(REAL(nx_ - 1, WP), 1.0_WP)
    dy = (ylim_(2) - ylim_(1)) / MAX(REAL(ny_ - 1, WP), 1.0_WP)

    ! Fill coordinate arrays
    do j = 1, nx_
        xs(j) = xlim_(1) + REAL(j - 1, WP) * dx
    end do
    do i = 1, ny_
        ys(i) = ylim_(1) + REAL(i - 1, WP) * dy
    end do

    ! Build observers (meshgrid equivalent)
    idx = 1
    do i = 1, ny_
        y_val = ys(i)
        do j = 1, nx_
            x_val = xs(j)
            observers_(idx, 1) = x_val
            observers_(idx, 2) = y_val
            observers_(idx, 3) = z_
            idx = idx + 1
        end do
    end do

    CALL self % check_observers_distances(observers_, 0.1_WP)

    write(*, '(A)') ''
    write(*, '(A, I0, A, I0, A, F0.2, A, F0.2, A, F0.2, A)') &
        'Computing XY slice (nx=', nx_, ', ny=', ny_, ', z=', z_, &
        ' m, center=(', center_(1), ', ', center_(2), ') m)...'

    CALL self % compute_pressure(observers_, block_size_, p)
    if (self % debug) write(*, '(A, F0.2)') 'SliceXY Pressure Norm: ', norm2(abs(p))

    ! Reshape pressure and observers to 3D (nx, ny, 3) and (Nfreqs, nx, ny)
    ALLOCATE(p_reshaped(SIZE(p,1), nx_, ny_))
    ALLOCATE(obs_reshaped(nx_, ny_, 3))

    idx = 1
    do i = 1, ny_
        do j = 1, nx_
            obs_reshaped(j, i, :) = observers_(idx, :)
            p_reshaped(:, j, i)   = p(:, idx)
            idx = idx + 1
        end do
    end do

    ! Save results
    CALL self % save_acoustics("Freqs", self % turbines(1) % Freqs)
    CALL self % save_acoustics("P_slicexy", p_reshaped)
    CALL self % save_acoustics("Obs_slicexy", obs_reshaped)
    CALL self % save_acoustics("X_slicexy", xs)
    CALL self % save_acoustics("Y_slicexy", ys)
    CALL self % save_acoustics("Z_slicexy", z_)
    CALL self % save_acoustics("Center_slicexy", center_)
    CALL self % save_acoustics("Nx_slicexy", nx_)
    CALL self % save_acoustics("Ny_slicexy", ny_)

    print*, ''
    print*, '   ---> SliceXY data saved at ', self % save_path

END SUBROUTINE run_sliceXY


SUBROUTINE run_sliceXZ(self, y, nx, nz, xlim, zlim, block_size)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: y
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nx
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nz
    REAL(WP)    , INTENT(IN), OPTIONAL  :: xlim(2)
    REAL(WP)    , INTENT(IN), OPTIONAL  :: zlim(2)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! Local variables
    COMPLEX(WP), ALLOCATABLE :: p(:,:)
    COMPLEX(WP), ALLOCATABLE :: p_reshaped(:,:,:)   ! (Nfreqs, nx, nz)
    REAL(WP)   , ALLOCATABLE :: observers_(:,:)
    REAL(WP)   , ALLOCATABLE :: obs_reshaped(:,:,:) ! (nx, nz, 3)
    REAL(WP)   , ALLOCATABLE :: xs(:), zs(:)        ! Coordinate arrays
    REAL(WP)                 :: y_, xlim_(2), zlim_(2), dx, dz, x_val, z_val
    INTEGER(I32)             :: nx_, nz_, block_size_, i, j, idx

    ! Defaults
    y_ = self % turbines(1) % AxisPos(2)         ; if (PRESENT(y))          y_          = y
    nx_ = 26_I32                                 ; if (PRESENT(nx))         nx_         = nx
    nz_ = 26_I32                                 ; if (PRESENT(nz))         nz_         = nz
    xlim_ = [-500.0_WP, 500.0_WP]                ; if (PRESENT(xlim))       xlim_       = xlim
    zlim_ = [0.0_WP, -self % turbines(1) % Depth]; if (PRESENT(zlim))       zlim_       = zlim
    block_size_ = 128_I32                        ; if (PRESENT(block_size)) block_size_ = block_size

    ALLOCATE(observers_(nx_ * nz_, 3))
    ALLOCATE(xs(nx_), zs(nz_))
    dx = (xlim_(2) - xlim_(1)) / MAX(REAL(nx_ - 1, WP), 1.0_WP)
    dz = (zlim_(2) - zlim_(1)) / MAX(REAL(nz_ - 1, WP), 1.0_WP)

    ! Fill coordinate arrays
    do j = 1, nx_
        xs(j) = xlim_(1) + REAL(j - 1, WP) * dx
    end do
    do i = 1, nz_
        zs(i) = zlim_(1) + REAL(i - 1, WP) * dz
    end do

    ! Build observers
    idx = 1
    do i = 1, nz_
        z_val = zs(i)
        do j = 1, nx_
            x_val = xs(j)
            observers_(idx, 1) = x_val
            observers_(idx, 2) = y_
            observers_(idx, 3) = z_val
            idx = idx + 1
        end do
    end do

    CALL self % check_observers_distances(observers_, 0.1_WP)

    write(*, '(A)') ''
    write(*, '(A, I0, A, I0, A, F0.2, A, F0.2, A, F0.2, A)') &
        'Computing XZ slice (nx=', nx_, ', nz=', nz_, ', y=', y_, &
        ' m, z from ', zlim_(1), ' to ', zlim_(2), ' m)...'

    CALL self % compute_pressure(observers_, block_size_, p)
    if (self % debug) write(*, '(A, F0.2)') 'SliceXZ Pressure Norm: ', norm2(abs(p))

    ! Reshape pressure and observers to 3D (nx, nz, 3) and (Nfreqs, nx, nz)
    ALLOCATE(p_reshaped(SIZE(p,1), nx_, nz_))
    ALLOCATE(obs_reshaped(nx_, nz_, 3))

    idx = 1
    do i = 1, nz_
        do j = 1, nx_
            obs_reshaped(j, i, :) = observers_(idx, :)
            p_reshaped(:, j, i)   = p(:, idx)
            idx = idx + 1
        end do
    end do

    ! Save results
    CALL self % save_acoustics("Freqs", self % turbines(1) % Freqs)
    CALL self % save_acoustics("P_slicexz", p_reshaped)
    CALL self % save_acoustics("Obs_slicexz", obs_reshaped)
    CALL self % save_acoustics("X_slicexz", xs)
    CALL self % save_acoustics("Z_slicexz", zs)
    CALL self % save_acoustics("Y_slicexz", y_)
    CALL self % save_acoustics("Nx_slicexz", nx_)
    CALL self % save_acoustics("Nz_slicexz", nz_)

    print*, ''
    print*, '   ---> SliceXZ data saved at ', self % save_path

END SUBROUTINE run_sliceXZ


SUBROUTINE dummy_run_sphere(self, r, n_theta, nz, center, block_size)
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
    REAL(WP)    , INTENT(IN), OPTIONAL  :: r
    INTEGER(I32), INTENT(IN), OPTIONAL  :: n_theta
    INTEGER(I32), INTENT(IN), OPTIONAL  :: nz
    REAL(WP)    , INTENT(IN), OPTIONAL  :: center(3)
    INTEGER(I32), INTENT(IN), OPTIONAL  :: block_size

    ! local variables not to have interpreter warnings
    REAL(WP) :: r_, center_(3)
    INTEGER(I32) ::n_theta_, nz_, block_size_

    
    r_ = 30.0_WP;          if (PRESENT(r))          r_ = r
    n_theta_ = 72_I32;     if (PRESENT(n_theta))    n_theta_ = n_theta
    nz_ = 20_I32;          if (PRESENT(nz))         nz_ = nz
    block_size_ = 128_I32; if (PRESENT(block_size)) block_size_ = block_size
    center_ = 0.0_WP;      if (PRESENT(center))     center_ = center

    ! Does nothing only skips
    print*, '  [Skipped] Spherical radiation not available for: ', trim(self % solver_name)

END SUBROUTINE dummy_run_sphere


SUBROUTINE run_all(self, &
    spectrums_observers, spectrums_z_obs, spectrums_block_size, &
    polar_r, polar_z, polar_n_theta, polar_center, polar_block_size, &
    cylinder_r, cylinder_n_theta, cylinder_nz, cylinder_center, cylinder_block_size, &
    decay_distance, decay_n_points, decay_z, decay_logspace, decay_block_size, &
    line_p1, line_p2, line_n_points, line_logspace, line_block_size, &
    sliceXY_z, sliceXY_nx, sliceXY_ny, sliceXY_xlim, sliceXY_ylim, sliceXY_center, sliceXY_block_size, &
    sliceXZ_y, sliceXZ_nx, sliceXZ_nz, sliceXZ_xlim, sliceXZ_zlim, sliceXZ_block_size, &
    sphere_r, sphere_n_theta, sphere_nz, sphere_center, sphere_block_size)
        
    CLASS(AcousticSolver_t), INTENT(INOUT) :: self
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
    CALL self % run_sphere(sphere_r, sphere_n_theta, sphere_nz, sphere_center, sphere_block_size)

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


END MODULE AcousticSolver