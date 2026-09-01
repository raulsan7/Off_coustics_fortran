MODULE AnalyticalNormalModes

USE omp_lib
USE IOUtils, ONLY: save_to_hdf5
USE Kinds, ONLY: WP, I32, PI, I1
USE WindTurbine, ONLY: WindTurbine_t
USE AcousticSolver, ONLY: AcousticSolver_t, save_turbine_params

IMPLICIT NONE

PRIVATE
PUBLIC :: AnalyticalNormalModes_t


! ------------------------
! Main Type Definition
! ------------------------
TYPE, EXTENDS(AcousticSolver_t) :: AnalyticalNormalModes_t

    INTEGER(I32) :: Nmodes    = 0           ! [-] Number of normal modes to retain
    REAL(WP)     :: H         = 30.0_WP     ! [m] Height of the waveguide (distance between the two boundaries)
    REAL(WP)     :: Upper_HBC = 0.0_WP      ! [m] Upper boundary height for the HBC
    REAL(WP)     :: Lower_HBC = -30.0_WP    ! [m] Lower boundary height for the HBC
    LOGICAL      :: verbose   = .false.     ! [-] Verbose output flag

    CONTAINS
    ! --- Public type-bound procedures --- !
    GENERIC   :: init        => init_array, init_single
    PROCEDURE :: init_array  => init_AnalyticalNormalModes
    PROCEDURE :: init_single => init_AnalyticalNormalModes_single
    ! PROCEDURE :: print_mode_summary
    ! PROCEDURE :: Psi_m
    ! PROCEDURE :: dPsi_m_dz
    ! PROCEDURE :: dipolar_pressure_NM

    ! --- Overridden deferred procedures --- !
    PROCEDURE :: get_name         => get_name_AnalyticalNormalModes
    PROCEDURE :: compute_pressure => compute_pressure_AnalyticalNormalModes
    PROCEDURE :: save_parameters  => save_parameters_AnalyticalNormalModes
    ! PROCEDURE :: save_acoustics   => save_acoustics_AnalyticalNormalModes

END TYPE AnalyticalNormalModes_t


CONTAINS

! ---------- DEFERRED ---------- !
SUBROUTINE init_AnalyticalNormalModes(self, turbines, Nmodes, c_wat, rho_wat, &
                                          Upper_HBC, Lower_HBC, eps, p_ref, &
                                          cluster, verbose, debug, name)
    
    CLASS(AnalyticalNormalModes_t), INTENT(INOUT) :: self
    CLASS(WindTurbine_t), INTENT(IN), OPTIONAL    :: turbines(:)
    INTEGER(I32)        , INTENT(IN), OPTIONAL    :: Nmodes
    REAL(WP)            , INTENT(IN), OPTIONAL    :: c_wat, rho_wat
    REAL(WP)            , INTENT(IN), OPTIONAL    :: Upper_HBC, Lower_HBC
    REAL(WP)            , INTENT(IN), OPTIONAL    :: eps, p_ref
    LOGICAL             , INTENT(IN), OPTIONAL    :: cluster, verbose, debug
    CHARACTER(len=*)    , INTENT(IN), OPTIONAL    :: name

    ! Local variables
    CHARACTER(len=512) :: name_

    ! Defaults
    if (PRESENT(Nmodes))    self % Nmodes    = Nmodes
    if (PRESENT(c_wat))     self % c         = c_wat
    if (PRESENT(rho_wat))   self % rho       = rho_wat
    if (PRESENT(Upper_HBC)) self % Upper_HBC = Upper_HBC
    if (PRESENT(Lower_HBC)) self % Lower_HBC = Lower_HBC
    if (PRESENT(eps))       self % eps       = eps
    if (PRESENT(p_ref))     self % p_ref     = p_ref
    if (PRESENT(cluster))   self % cluster   = cluster
    if (PRESENT(verbose))   self % verbose   = verbose
    if (PRESENT(debug))     self % debug     = debug

    ! Build save path
    self % save_path = trim(self % turbines(1) % save_path) // trim(name) // ".hdf5"
    self % H = abs(self % Upper_HBC - self % Lower_HBC)

    if (self % Upper_HBC <= self % Lower_HBC) error stop "MethodImages % init: Upper_HBC must be strictly greater than Lower_HBC"
    
    if (PRESENT(turbines)) then
        if (ALLOCATED(self % turbines)) DEALLOCATE(self % turbines)
        ALLOCATE(self % turbines, source = turbines)
    end if

    print*, ''
    print*, 'Selected method: ', self % get_name()
    ! if (self % verbose) CALL self % print_mode_summary()
    ! CALL self % save_parameters()

END SUBROUTINE init_AnalyticalNormalModes


SUBROUTINE init_AnalyticalNormalModes_single(self, turbine, Nmodes, c_wat, rho_wat, &
                                                 Upper_HBC, Lower_HBC, eps, p_ref, &
                                                 cluster, verbose, debug, name)
        CLASS(AnalyticalNormalModes_t), INTENT(INOUT) :: self
        CLASS(WindTurbine_t) , INTENT(IN)            :: turbine
        INTEGER(I32)         , INTENT(IN), OPTIONAL  :: Nmodes
        REAL(WP)             , INTENT(IN), OPTIONAL  :: c_wat, rho_wat
        REAL(WP)             , INTENT(IN), OPTIONAL  :: Upper_HBC, Lower_HBC
        REAL(WP)             , INTENT(IN), OPTIONAL  :: eps, p_ref
        LOGICAL              , INTENT(IN), OPTIONAL  :: cluster, verbose, debug
        CHARACTER(len=*)     , INTENT(IN)            :: name

        if (ALLOCATED(self % turbines)) DEALLOCATE(self % turbines)
        ALLOCATE(self % turbines(1), source = turbine)

        CALL init_AnalyticalNormalModes(self, Nmodes=Nmodes, c_wat=c_wat, rho_wat=rho_wat, &
                                        Upper_HBC=Upper_HBC, Lower_HBC=Lower_HBC, eps=eps, p_ref=p_ref, &
                                        cluster=cluster, verbose=verbose, debug=debug, name=name)

END SUBROUTINE init_AnalyticalNormalModes_single


FUNCTION get_name_AnalyticalNormalModes(self) RESULT(name)
    CLASS(AnalyticalNormalModes_t), INTENT(INOUT) :: self
    CHARACTER(len=:)     , ALLOCATABLE   :: name
    name = "Analytical Normal Modes"; self % solver_name = trim(name)
END FUNCTION get_name_AnalyticalNormalModes


SUBROUTINE compute_pressure_AnalyticalNormalModes(self, observers, block_size, total_pressure)
    CLASS(AnalyticalNormalModes_t) , INTENT(INOUT) :: self
    REAL(WP) , INTENT(IN) :: observers(:,:) 
    INTEGER(I32) , INTENT(IN), OPTIONAL :: block_size 
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT) :: total_pressure(:,:) 

    INTEGER(I32) :: nobs, nf, nturb, nt, bs, n_chunks, ichunk
    INTEGER(I32) :: istart, iend, nb, i
    REAL(WP), ALLOCATABLE :: freqs(:)
    COMPLEX(WP), ALLOCATABLE :: p_block(:,:)
    COMPLEX(WP), ALLOCATABLE :: F_corr(:,:,:)

    bs = 256_I32
    if (PRESENT(block_size)) bs = block_size

    ! nobs = SIZE(observers, 1)
    ! if (nobs <= 0) then
    !     allocate(total_pressure(0,0))
    !     return
    ! end if

    ! if (SIZE(observers, 2) /= 3) error stop "compute_pressure_AnalyticalNormalModes: observers shape must be (N,3)"

    ! nturb = SIZE(self % turbines)
    ! freqs = self % turbines(1) % Freqs
    ! nf = SIZE(freqs)
    ! allocate(total_pressure(nf, nobs))
    ! total_pressure = (0.0_WP, 0.0_WP)

    ! if (.not. self % cluster) bs = 1

    ! if (self % cluster) then
    !     do nt = 1, nturb
    !         if (nturb > 1) then
    !             print '(A,I0,A,I0,A)', 'Turbine ', nt, '/', nturb, ':'
    !             flush(6)
    !         end if
    !         ! Aplicamos corrección por impedancia como en MethodImages
    !         F_corr = self % turbines(nt) % get_impedance_corrected_force(self % c)
    !         n_chunks = (nobs + bs - 1) / bs

    !         !$OMP PARALLEL DO DEFAULT(SHARED) &
    !         !$OMP PRIVATE(ichunk, istart, iend, nb, p_block)&
    !         !$OMP SHARED(freqs, n_chunks, nobs, total_pressure, self, nt, F_corr, observers)
    !         do ichunk = 0, n_chunks - 1
    !             istart = ichunk * bs + 1
    !             iend = min(istart + bs - 1, nobs)
    !             nb = iend - istart + 1

    !             call self % dipolar_pressure_NM(observers(istart:iend, :), freqs, &
    !                                             self % turbines(nt) % x, F_corr, p_block)

    !             total_pressure(:, istart:iend) = total_pressure(:, istart:iend) + p_block(:, 1:nb)
    !             if (ALLOCATED(p_block)) DEALLOCATE(p_block)
    !         end do
    !         !$OMP END PARALLEL DO
    !     end do
    ! else
    !     do nt = 1, nturb
    !         if (nturb > 1) then
    !             print '(A,I0,A,I0,A)', 'Turbine ', nt, '/', nturb, ':'
    !         end if
    !         F_corr = self % turbines(nt) % get_impedance_corrected_force(self % c)
    !         do i = 1, nobs
    !             call self % dipolar_pressure_NM(observers(i:i, :), freqs, &
    !                                             self % turbines(nt) % x, F_corr, p_block)
    !             total_pressure(:, i) = total_pressure(:, i) + p_block(:, 1)
    !             if (ALLOCATED(p_block)) DEALLOCATE(p_block)

    !             if (mod(i, 100) == 0) then
    !                 print '(A,I0,A,I0)', '    Progress: ', i, '/', nobs
    !                 flush(6)
    !             end if
    !         end do
    !     end do
    ! end if

END SUBROUTINE compute_pressure_AnalyticalNormalModes


SUBROUTINE save_parameters_AnalyticalNormalModes(self)
    CLASS(AnalyticalNormalModes_t), INTENT(IN) :: self

    ! Local variables
    CHARACTER(len=512) :: save_path
    
    save_path = trim(self % save_path)

    ! Solver
    CALL save_to_hdf5(save_path, "Nmodes", self % Nmodes)
    CALL save_to_hdf5(save_path, "H", self % H)
    CALL save_to_hdf5(save_path, "c_wat"    , self % c)
    CALL save_to_hdf5(save_path, "rho_wat"  , self % rho)
    CALL save_to_hdf5(save_path, "Upper_HBC", self % Upper_HBC)
    CALL save_to_hdf5(save_path, "Lower_HBC", self % Lower_HBC)
    CALL save_to_hdf5(save_path, "p_ref"    , self % p_ref)
    CALL save_to_hdf5(save_path, "Nturb"    , size(self % turbines))

    CALL save_turbine_params(self)

END SUBROUTINE save_parameters_AnalyticalNormalModes


END MODULE AnalyticalNormalModes