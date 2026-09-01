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
    PROCEDURE :: print_mode_summary
    PROCEDURE :: Psi_kzm
    PROCEDURE :: dPsi_kzm_dz
    PROCEDURE :: dipolar_pressure_NM

    ! --- Overridden deferred procedures --- !
    PROCEDURE :: get_name         => get_name_AnalyticalNormalModes
    PROCEDURE :: compute_pressure => compute_pressure_AnalyticalNormalModes
    PROCEDURE :: save_parameters  => save_parameters_AnalyticalNormalModes

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
    REAL(WP)           :: f_max
    INTEGER(I32)       :: m_prop, near_field_buffer

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

    if (.not. ALLOCATED(self % turbines)) then
        error stop "AnalyticalNormalModes % init: No turbines defined."
    end if

    ! Compute default number of modes if not specified
    if ( self % Nmodes <= 0) then
        f_max = maxval(self % turbines(1) % Freqs)
        m_prop = int(floor(0.5_WP * (4.0_WP * self % H * f_max / self % c + 1.0_WP)))
        near_field_buffer = 5
        self % Nmodes = max(1_I32, m_prop + near_field_buffer)
    end if

    print*, ''
    print*, 'Selected method: ', self % get_name()
    if (self % verbose) CALL self % print_mode_summary()
    CALL self % save_parameters()

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
    CHARACTER(len=:)              , ALLOCATABLE   :: name
    name = "Analytical Normal Modes"; self % solver_name = trim(name)
END FUNCTION get_name_AnalyticalNormalModes


SUBROUTINE compute_pressure_AnalyticalNormalModes(self, observers, block_size, total_pressure)
    CLASS(AnalyticalNormalModes_t) , INTENT(INOUT)        :: self
    REAL(WP)                       , INTENT(IN)           :: observers(:,:) 
    INTEGER(I32)                   , INTENT(IN), OPTIONAL :: block_size 
    COMPLEX(WP), ALLOCATABLE       , INTENT(OUT)          :: total_pressure(:,:) 

    INTEGER(I32) :: nobs, nf, nturb, nt, bs, n_chunks, ichunk
    INTEGER(I32) :: istart, iend, nb, i
    REAL(WP), ALLOCATABLE :: freqs(:)
    COMPLEX(WP), ALLOCATABLE :: p_block(:,:)
    COMPLEX(WP), ALLOCATABLE :: F_corr(:,:,:)

    bs = 256_I32
    if (PRESENT(block_size)) bs = block_size

    nobs = SIZE(observers, 1)
    if (nobs <= 0) then
        allocate(total_pressure(0,0))
        return
    end if

    if (SIZE(observers, 2) /= 3) error stop "compute_pressure_AnalyticalNormalModes: observers shape must be (N,3)"

    nturb = SIZE(self % turbines)
    freqs = self % turbines(1) % Freqs
    nf = SIZE(freqs)

    allocate(total_pressure(nf, nobs))
    total_pressure = (0.0_WP, 0.0_WP)

    if (.not. self % cluster) bs = 1

    if (self % cluster) then
        do nt = 1, nturb

            if (nturb > 1) then
                print '(A,I0,A,I0,A)', 'Turbine ', nt, '/', nturb, ':'
                flush(6)
            end if
            
            n_chunks = (nobs + bs - 1) / bs

            !$OMP PARALLEL DO DEFAULT(SHARED) &
            !$OMP PRIVATE(ichunk, istart, iend, nb, p_block)&
            !$OMP SHARED(freqs, n_chunks, nobs, total_pressure, self, nt, F_corr, observers)
            do ichunk = 0, n_chunks - 1
                istart = ichunk * bs + 1
                iend   = min(istart + bs - 1, nobs)
                nb     = iend - istart + 1
                
                ! Call kernel for this chunk (p_block will be allocated by the routine)
                call self % dipolar_pressure_NM(observers(istart:iend, :), freqs, &
                                                self % turbines(nt) % x, self % turbines(nt) % F, p_block)
                
                ! Accumulate results into the global array (non-overlapping slices -> no race)
                total_pressure(:, istart:iend) = total_pressure(:, istart:iend) + p_block(:, 1:nb)

                if (ALLOCATED(p_block)) DEALLOCATE(p_block)
            end do
            !$OMP END PARALLEL DO
        end do
    else
        do nt = 1, nturb

            if (nturb > 1) then
                print '(A,I0,A,I0,A)', 'Turbine ', nt, '/', nturb, ':'
            end if

            ! Serial: loop observers (keep 2D slice for compatibility)
            do i = 1, nobs
                call self % dipolar_pressure_NM(observers(i:i, :), freqs, &
                                                self % turbines(nt) % x, self % turbines(nt) % F, p_block)
                total_pressure(:, i) = total_pressure(:, i) + p_block(:, 1)

                if (ALLOCATED(p_block)) DEALLOCATE(p_block)

                ! Print progress every 100 observers
                if (mod(i, 100) == 0) then
                    print '(A,I0,A,I0)', '    Progress: ', i, '/', nobs
                    flush(6)
                end if
            end do
        end do
    end if

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


! ---------- HELPERS ---------- !
SUBROUTINE print_mode_summary(self)
    CLASS(AnalyticalNormalModes_t), INTENT(IN) :: self

    ! Local variables
    REAL(WP)          :: fmin, fmax, k_medium, lambda_zm, fc, kr_sq, k_zm_val, krm_val
    COMPLEX(WP)       :: kzm, krm
    INTEGER(I32)      :: i, m_idx, line_width
    CHARACTER(len=20) :: lambda_str, status
    INTEGER(I32), ALLOCATABLE :: modes_to_show(:)


    fmin = minval(self % turbines(1) % Freqs)
    fmax = maxval(self % turbines(1) % Freqs)
    k_medium = 2.0_WP * PI * fmax / self % c
    line_width = 115

    print*, repeat('=', line_width)
    write(*, '(A, T43, A)') " ", "ANALYTICAL NORMAL MODES SETUP"
    print *, repeat("=", line_width)
    write(*, '(A, F10.2, A)') "  Water Depth (H):      ", self % H, " m"
    write(*, '(A, F10.2, A)') "  Speed of Sound (c):   ", self % c, " m/s"
    write(*, '(A, F10.2, A, F10.2, A)') "  Frequency Band:       ", fmin, " Hz - ", fmax, " Hz"
    write(*, '(A, I5)') "  Retained Modes (m):   ", self % Nmodes
    write(*, '(A, F10.2, A)') "  * Note: k_rm and Horiz. Wavelength are evaluated at f_max (", fmax, " Hz)"
    print *, repeat("-", line_width)
    write(*, '(A5, A2, A18, A3, A12, A3, A13, A3, A20, A3, A13, A3, A11)') &
        "Mode", " |", "k_zm [rad/m]", " |", "Vert. WL [m]", " |", "f_cutoff [Hz]", " |", &
        "k_rm [rad/m]", " |", "Horiz. WL [m]", " |", "Status"
    print *, repeat("-", line_width)

    if (self % Nmodes <= 20) then
        allocate(modes_to_show(self % Nmodes))
        do i = 1, self % Nmodes
            modes_to_show(i) = i
        end do
    else
        allocate(modes_to_show(16))
        do i = 1, 10
            modes_to_show(i) = i
        end do
        modes_to_show(11) = -1 
        do i = 12, 16
            modes_to_show(i) = self % Nmodes - (16 - i)
        end do
    end if

    do i = 1, size(modes_to_show)
        m_idx = modes_to_show(i)
        if (m_idx == -1) then
            write(*, '(A5, A2, A18, A3, A12, A3, A13, A3, A20, A3, A13, A3, A11)') &
                "...", " |", "...", " |", "...", " |", "...", " |", "...", " |", "...", " |", "..."
            cycle
        end if

        k_zm_val = (2.0_WP * REAL(m_idx, WP) - 1.0_WP) * PI / (2.0_WP * self % H)
        kzm = cmplx(k_zm_val, 0.0_WP, kind=WP)
        lambda_zm = 4.0_WP * self % H / (2.0_WP * REAL(m_idx, WP) - 1.0_WP)
        fc = (2.0_WP * REAL(m_idx, WP) - 1.0_WP) * self % c / (4.0_WP * self % H)

        kr_sq = k_medium**2 - k_zm_val**2
        if (kr_sq >= 0.0_WP) then
            krm_val = sqrt(kr_sq)
            krm = cmplx(krm_val, 0.0_WP, kind=WP)
            write(lambda_str, '(F10.2)') 2.0_WP * PI / krm_val
            status = "Propagating"
        else
            krm_val = sqrt(abs(kr_sq))
            krm = cmplx(0.0_WP, krm_val, kind=WP)
            lambda_str = "N/A (Decays)"
            status = "Evanescent"
        end if

        write(*, '(I5, A2, F12.4, A6, A2, F12.2, A2, F13.2, A2, F12.4, A1, F12.4, A1, A2, A13, A2, A11)') &
            m_idx, " |", real(kzm, WP), "+0.0000j", " |", lambda_zm, " |", fc, " |", &
            real(krm, WP), "+", aimag(krm), "j", " |", trim(adjustl(lambda_str)), " |", status
    end do
    print *, repeat("=", line_width)
    print *, ""




END SUBROUTINE print_mode_summary


! ---------- METHOD ---------- !
FUNCTION Psi_kzm(self, z, m, kzm) RESULT(psi)
    CLASS(ANalyticalNormalModes_t), INTENT(IN) :: self
    REAL(WP)    , INTENT(IN)           :: z(:)
    INTEGER(I32), INTENT(IN), OPTIONAL :: m
    REAL(WP)    , INTENT(IN), OPTIONAL :: kzm

    ! Local variables
    REAL(WP), ALLOCATABLE :: psi(:)

    ALLOCATE(psi(size(z)))

    if (PRESENT(m) .and. .not. PRESENT(kzm)) then
        psi = sqrt(2.0_WP * self % rho / self % H) * &
              sin((2.0_WP * REAL(m, WP) - 1.0_WP) * PI * z / (2.0_WP * self % H))
    else if (PRESENT(kzm)) then
        psi = sqrt(2.0_WP * self % rho / self % H) * sin(kzm * z)
    else
        error stop "Psi_kzm: Either m or kzm must be provided."
    end if
    
END FUNCTION Psi_kzm


FUNCTION dPsi_kzm_dz(self, z, m, kzm) RESULT(dpsi)
    CLASS(AnalyticalNormalModes_t), INTENT(IN) :: self
    REAL(WP)    , INTENT(IN)           :: z(:)
    INTEGER(I32), INTENT(IN), OPTIONAL :: m
    REAL(WP)    , INTENT(IN), OPTIONAL :: kzm

    ! Local variables
    REAL(WP), ALLOCATABLE :: dpsi(:)
    REAL(WP)              :: k_zm
    
    if (PRESENT(kzm)) then
        k_zm = kzm
    else if (PRESENT(m)) then
        k_zm = (2.0_WP * REAL(m, WP) - 1.0_WP) * PI / (2.0_WP * self % H)
    else
        error stop "dPsi_kzm_dz: Either m or kzm must be provided."
    end if

    allocate(dpsi(size(z)))
    dpsi = sqrt(2.0_WP * self % rho / self % H) * k_zm * cos(k_zm * z)

END FUNCTION dPsi_kzm_dz


SUBROUTINE dipolar_pressure_NM(self, obs_pos, freqs, nodes_pos, force, p_out)
    USE MathUtils, ONLY: compute_hankel_complex
    CLASS(AnalyticalNormalModes_t), INTENT(IN) :: self
    REAL(WP)                      , INTENT(IN) :: obs_pos(:,:)
    REAL(WP)                      , INTENT(IN) :: freqs(:)
    REAL(WP)                      , INTENT(IN) :: nodes_pos(:,:)
    COMPLEX(WP)                   , INTENT(IN) :: force(:,:,:)      ! (Nfreqs, Nnodes, 3)
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT)      :: p_out(:,:) 

    ! Local variables
    INTEGER(I32) :: Nfreqs, Nnodes, Nob, j, f, k, m
    REAL(WP)     :: epsv, k_zm, Y, Z_real, kappa
    COMPLEX(WP)  :: k_rm, H0, H1, term_z, dx_part, dy_part, term_xy, p_m
    COMPLEX(WP)  :: fx, fy, fz
    REAL(WP), ALLOCATABLE :: dx(:,:), dy(:,:), r(:,:)
    REAL(WP), ALLOCATABLE :: k_wave(:)
    REAL(WP), ALLOCATABLE :: z_obs(:), zs(:)
    REAL(WP), ALLOCATABLE :: psi_z(:), psi_zs(:), dpsi_zs(:)

    Nfreqs = size(freqs)
    Nnodes = size(nodes_pos, 1)
    Nob = size(obs_pos, 1)
    epsv = self % eps

    ALLOCATE(p_out(Nfreqs, Nob), k_wave(Nfreqs))
    p_out = (0.0_WP, 0.0_WP)
    k_wave = 2.0_WP * PI * freqs / self % c

    ALLOCATE(dx(Nob, Nnodes), dy(Nob, Nnodes), r(Nob, Nnodes))
    ALLOCATE(z_obs(Nob), zs(Nnodes))
    z_obs = obs_pos(:, 3)
    zs    = nodes_pos(:, 3)

    ! Precompute horizontal distance differences
    do j = 1, Nnodes
        do k = 1, Nob
            dx(k, j) = obs_pos(k, 1) - nodes_pos(j, 1)
            dy(k, j) = obs_pos(k, 2) - nodes_pos(j, 2)
            r(k, j)  = sqrt(dx(k, j)**2 + dy(k, j)**2)
        end do
    end do

    ALLOCATE(psi_z(Nob), psi_zs(Nnodes), dpsi_zs(Nnodes))

    ! Loop over frequencies
    do f = 1, Nfreqs
        ! Normal modes superposition
        do m = 1, self % Nmodes
            k_zm = (2.0_WP * REAL(m, WP) - 1.0_WP) * PI / (2.0_WP * self % H)

            ! Compute krm avoiding external libraries
            if (k_wave(f) >= k_zm) then
                k_rm = cmplx(sqrt(k_wave(f)**2 - k_zm**2), 0.0_WP, kind=WP)
            else
                kappa = sqrt(k_zm**2 - k_wave(f)**2)
                k_rm = cmplx(0.0_WP, kappa, kind=WP)
            end if

            ! Evaluate normal modes
            psi_z   = self % Psi_kzm    (z_obs, kzm=k_zm)
            psi_zs  = self % Psi_kzm    (zs   , kzm=k_zm)
            dpsi_zs = self % dPsi_kzm_dz(zs   , kzm=k_zm)

            ! Loop over nodes and receivers
            do j = 1, Nnodes
                fx = force(f, j, 1)
                fy = force(f, j, 2)
                fz = force(f, j, 3)

                do k = 1, Nob
                    ! Compute hankel functions
                    CALL compute_hankel_complex(k_rm, r(k,j), H0, H1)

                    ! Vertical component
                    term_z = fz * (dpsi_zs(j) / self % rho) * H0

                    ! Horizontal components
                    dx_part = fx * (dx(k,j)/r(k,j))
                    dy_part = fy * (dy(k,j)/r(k,j))
                    term_xy = (psi_zs(j) * k_rm / self % rho) * H1 * (dx_part + dy_part)

                    ! Linear superpositions of node m on receiver k
                    p_m = psi_z(k) * (term_z + term_xy)
                    p_out(f, k) = p_out(f, k) + p_m
                end do
            end do
        end do
    end do

    ! Global coefficient
    p_out = -0.25_WP * I1 * p_out


END SUBROUTINE dipolar_pressure_NM


END MODULE AnalyticalNormalModes