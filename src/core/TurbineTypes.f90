MODULE TurbineTypes

USE Kinds, ONLY: WP, I32, PI
USE WindTurbine

IMPLICIT NONE 
PRIVATE

PUBLIC:: DTU10MWMonopile

TYPE, EXTENDS(WindTurbine_t) :: DTU10MWMonopile
    ! Physical and structural parameters specific to the DTU 10 MW Monopile
    REAL(WP) :: D        = 9.0_WP      ! [m] Base monopile outer diameter
    REAL(WP) :: rho_wat  = 1025.0_WP   ! [kg/m^3] Seawater density
    REAL(WP) :: rho_mat  = 8500.0_WP   ! [kg/m^3] Material density (steel)
    REAL(WP) :: wet_area = 848.12_WP   ! [m^2] Wetted area

    ! Path to rotor speed curve CSV file
    CHARACTER(len=:), ALLOCATABLE :: path_rpm

    CONTAINS
    PROCEDURE :: init                          => init_DTU10MWMonopile
    PROCEDURE :: compute_force                 => compute_force_DTU10MWMonopile
    PROCEDURE :: compute_mass                  => compute_mass_DTU10MWMonopile
    PROCEDURE :: filter_frequencies            => filter_frequencies_DTU10MWMonopile
    ! PROCEDURE :: get_impedance_corrected_force => get_impedance_corrected_force_DTU10MWMonopile


END TYPE DTU10MWMonopile


CONTAINS

! ---------- TURBINE MODEL 1 ---------- !
SUBROUTINE init_DTU10MWMonopile(self, debug, rootname, output_dir, save_dir, save_name, &
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

END SUBROUTINE init_DTU10MWMonopile


SUBROUTINE compute_force_DTU10MWMonopile(self, filter_freqs, verbose, skipf)
    USE MathUtils, ONLY: remove_duplicate_nodes, compute_rfft

    CLASS(DTU10MWMonopile), INTENT(INOUT)        :: self
    LOGICAL               , INTENT(IN), OPTIONAL :: filter_freqs
    LOGICAL               , INTENT(IN), OPTIONAL :: verbose
    INTEGER(I32)          , INTENT(IN), OPTIONAL :: skipf

    ! Local variables
    INTEGER(I32) :: Nm, Nn, Nt, Nnodes_wet, i, j, k, idx, Nfreqs
    REAL(WP)     :: dt, F_memory_mb
    LOGICAL     , ALLOCATABLE :: mask_duplicate(:), wet_nodes(:), mask_freqs_to_use(:)
    REAL(WP)    , ALLOCATABLE :: mass(:), added_mass(:), mass_effective(:)
    REAL(WP)    , ALLOCATABLE :: x_flat(:,:)         ! (Nm*Nn, 3)
    REAL(WP)    , ALLOCATABLE :: acc_flat(:,:,:)     ! (nt, Nm*Nn, 3)
    COMPLEX(WP) , ALLOCATABLE :: A(:,:,:)            ! FFT result (Nfreqs, Nnodes_wet, 3)
    REAL(WP)    , ALLOCATABLE :: freqs(:)
    INTEGER(I32), ALLOCATABLE :: idxs(:)
    LOGICAL      :: verbose_ = .false., filter_freqs_ = .true.
    INTEGER(I32) :: skipf_ = 1

    ! Defaults
    if (PRESENT(filter_freqs)) filter_freqs_ = filter_freqs
    if (PRESENT(verbose))      verbose_      = verbose
    if (PRESENT(skipf))        skipf_        = skipf

    ! Get dimensions
    Nm = self % Nmembers
    Nn = self % Nnodes
    nt = size(self % Time)
    dt = self % Time(2) - self % Time(1)

    ! Remove duplicate nodes shared between consecutive members
    mask_duplicate = remove_duplicate_nodes(self % x_all, delete_last=.true.)

    
    ! Reshape x_all and acc to 2D/3D preserving (member, node) ordering
    allocate(x_flat(Nm*Nn, 3))
    allocate(acc_flat(nt, Nm*Nn, 3))
    idx = 0
    do i = 1, Nm
        do j = 1, Nn
            idx = idx + 1
            x_flat(idx, :) = self % x_all(i, j, :)
        end do
    end do

    do i = 1, nt
        idx = 0
        do j = 1, Nm
            do k = 1, Nn
                idx = idx + 1
                acc_flat(i, idx, :) = self % acc(i, j, k, :)
            end do
        end do
    end do

    ! Apply duplicate mask: build integer index array from logical mask
    if (count(mask_duplicate) > 0) then
        idxs = pack([(i, i = 1, Nm*Nn)], mask_duplicate)
        x_flat = x_flat(idxs, :)
        acc_flat = acc_flat(:, idxs, :)
        DEALLOCATE(idxs)
    else ! No nodes selected: return zero-sized arrays
        x_flat = x_flat(0: -1, :)
        acc_flat = acc_flat(:, 0: -1, :)
    end if

    ! Remove dry nodes (z > 0): keep nodes with z <= 0
    wet_nodes = x_flat(:,3) <= 0.0_WP
    Nnodes_wet = count(wet_nodes)
    if (Nnodes_wet > 0) then
        idxs = pack([(i, i = 1, size(x_flat,1))], wet_nodes)
        x_flat = x_flat(idxs, :)
        acc_flat = acc_flat(:, idxs, :)
        DEALLOCATE(idxs)
    else
        x_flat = x_flat(0: -1, :)
        acc_flat = acc_flat(:, 0: -1, :)
    end if

    ! Store in a new variable
    self % x = x_flat

    ! Compute mass properties
    CALL self % compute_mass(mass, added_mass)
    mass_effective = mass + added_mass
    DEALLOCATE(mass, added_mass)

    ! Compute force via filter_freqs
    CALL compute_rfft(acc_flat, nt, dt, skipf=skipf_, remove_zero=.true., &
                      array_out=A, freqs=freqs)

    
    ! F = -A * mass (broadcast in frequency and direction)
    self % Freqs = freqs
    ALLOCATE(self % F(size(freqs), Nnodes_wet, 3))
    do i = 1, size(freqs)
        do j = 1, Nnodes_wet
            self % F(i, j, :) = - A(i, j, :) * mass_effective(j)
        end do
    end do

    ! Filter frequencies (optional)
    if (filter_freqs) then
        mask_freqs_to_use = self % filter_frequencies()
        idxs = pack([(i, i = 1, size(self%Freqs))], mask_freqs_to_use)

        self % Freqs = self % Freqs(idxs)
        self % F     = self % F(idxs, :, :)
    end if

    ! 6. Cleanup
    deallocate(self%acc, self%Time, A, freqs, mass_effective, mask_duplicate, wet_nodes)

    ! 7. Summary print (if verbose)
    Nfreqs = size(self%Freqs)
    F_memory_mb = real(storage_size(self%F) * size(self%F) / 8) / 1024.0_WP / 1024.0_WP
    if (verbose) then
        print '(A)', repeat('=', 70)
        print '(A)', 'FORCE COMPUTATION COMPLETED'
        print '(A)', repeat('-', 68)
        print '(A,I6)',   'Nodes (wet): ', Nnodes_wet
        print '(A,I6)',   'Frequencies: ', Nfreqs
        print '(A,F6.2,A)', 'Memory (F):  ', F_memory_mb, ' MB'
        print '(A)', repeat('=', 70)
    end if


END SUBROUTINE compute_force_DTU10MWMonopile


SUBROUTINE compute_mass_DTU10MWMonopile(self, mass, added_mass)
    USE MathUtils, ONLY: divide_span

    CLASS(DTU10MWMonopile), INTENT(INOUT) :: self
    REAL(WP), ALLOCATABLE, INTENT(OUT)    :: mass(:)
    REAL(WP), ALLOCATABLE, INTENT(OUT)    :: added_mass(:)

    ! Local variables
    REAL(WP) :: D, Depth, rho_mat, rho_wat
    INTEGER(I32) :: Nnodes, i
    REAL(WP), ALLOCATABLE :: z_nodes(:), ds(:), D_outer(:), D_inner(:)
    REAL(WP), ALLOCATABLE :: A_mass(:), A_added_mass(:)

    ! Retrieve parameters
    D       = self%D
    Depth   = self%Depth
    rho_mat = self%rho_mat
    rho_wat = self%rho_wat
    Nnodes  = size(self%x, 1)

    ! Generate vertical coordinate array linearly from -Depth to 0
    allocate(z_nodes(Nnodes))
    do i = 1, Nnodes
        z_nodes(i) = -Depth + real(i-1, WP) * Depth / real(Nnodes-1, WP)
    end do

    ! Nodal spacings
    ds = divide_span(z_nodes)

    ! Outer diameter constant
    allocate(D_outer(Nnodes))
    D_outer = D

    ! Inner diameter depending on elevation
    allocate(D_inner(Nnodes))
    do i = 1, Nnodes
        if (z_nodes(i) >= 4.0_WP .and. z_nodes(i) <= 10.0_WP) then
            D_inner(i) = 8.70_WP
        else if (z_nodes(i) >= -10.0_WP .and. z_nodes(i) < 4.0_WP) then
            D_inner(i) = 8.69_WP
        else if (z_nodes(i) < -10.0_WP) then
            D_inner(i) = 8.80_WP
        else
            D_inner(i) = 0.0_WP  ! should not happen
        end if
    end do

    ! Cross-sectional areas
    allocate(A_mass(Nnodes), A_added_mass(Nnodes))
    A_mass       = PI * ((D_outer/2.0_WP)**2 - (D_inner/2.0_WP)**2)
    A_added_mass = PI * (D_outer/2.0_WP)**2

    ! Distributed masses
    allocate(mass(Nnodes), added_mass(Nnodes))
    mass       = A_mass       * rho_mat * ds
    added_mass = A_added_mass * rho_wat * ds

    ! Cleanup local allocatables (optional, automatic on return)
    deallocate(z_nodes, ds, D_outer, D_inner, A_mass, A_added_mass)

END SUBROUTINE compute_mass_DTU10MWMonopile


FUNCTION filter_frequencies_DTU10MWMonopile(self) RESULT(mask)
    USE IOUtils, ONLY: read_curve
    USE MathUtils, ONLY: generate_timeseries_banded_sines, filter_non_usefull_freqs

    CLASS(DTU10MWMonopile), INTENT(INOUT) :: self
    LOGICAL, ALLOCATABLE                  :: mask(:)

    ! Local variables
    CHARACTER(len=20), ALLOCATABLE :: keys(:)
    REAL(WP)         , ALLOCATABLE :: freqs_amp(:,:)
    REAL(WP)         , ALLOCATABLE :: freqs_to_use(:)
    REAL(WP)                       :: rpm, freqs_over = 10.0_WP

    REAL(WP), ALLOCATABLE :: tmp_signal(:)

    ! Read rpm at current wind speed
    rpm = read_curve(self % path_rpm, self % WindSpeed)

    ! Generate drivetrain excitation spectrum
    CALL drivetrain10MW_excitation_spectrum(rpm=rpm, freqs_amp=freqs_amp, keys=keys)

    ! Get allowed frequencies
    CALL generate_timeseries_banded_sines(peaks=freqs_amp, keys=keys, t=self % Time, &
                                         zeta=0.02_WP, used_freqs=.true., a=tmp_signal, &
                                         freqs_out=freqs_to_use)

    ! Build mask
    print *, "DTU10MWMonopile.filter_frequencies(): FREQS_OVER IS HARDCODED TO 10.0 Hz"
    mask = filter_non_usefull_freqs(self % Freqs, freqs_to_use, freqs_over)

    DEALLOCATE(tmp_signal)

    DEALLOCATE(freqs_amp, keys, freqs_to_use)

END FUNCTION filter_frequencies_DTU10MWMonopile



! ---------- Drivetrain excitations ---------- !
SUBROUTINE drivetrain10MW_excitation_spectrum(rpm, damping, p_shaft, alpha_mesh, freqs_amp, keys)
    !> Constructs the physically-consistent excitation spectrum for the
    !> DTU 10 MW medium-speed drivetrain.
    REAL(WP), INTENT(IN), OPTIONAL :: rpm
    REAL(WP), INTENT(IN), OPTIONAL :: damping
    REAL(WP), INTENT(IN), OPTIONAL :: p_shaft
    REAL(WP), INTENT(IN), OPTIONAL :: alpha_mesh
    REAL(WP), ALLOCATABLE, INTENT(OUT) :: freqs_amp(:,:)        ! shape(N,2)
    CHARACTER(len=20), ALLOCATABLE, INTENT(OUT) :: keys(:)      ! shape(N)

    ! Local variables
    REAL(WP), PARAMETER :: min_rpm = 6.0_WP, max_rpm = 9.6_WP
    REAL(WP) :: factor, rpm_clipped, modal_sum

    ! Torsional eigenfrequencies (Hz)
    REAL(WP), PARAMETER :: torsional_modes(14) = [ &
        3.942_WP, 16.683_WP, 33.811_WP, 61.736_WP, 72.015_WP, &
        108.447_WP, 175.792_WP, 184.351_WP, 214.403_WP, &
        245.476_WP, 274.443_WP, 336.874_WP, 346.436_WP, 391.730_WP ]

    ! Shaft frequenciess at min and max rpm
    REAL(WP), PARAMETER :: shaft_min(13) = [ &
        0.1_WP, 0.2_WP, 0.3_WP, 0.6_WP, &
        0.322_WP, 0.644_WP, 0.966_WP, &
        1.383_WP, 2.766_WP, 4.149_WP, &
        5.017_WP, 10.034_WP, 15.051_WP ]
    REAL(WP), PARAMETER :: shaft_max(13) = [ &
        0.160_WP, 0.320_WP, 0.480_WP, 0.640_WP, &
        0.515_WP, 1.030_WP, 1.545_WP, &
        2.213_WP, 4.426_WP, 6.639_WP, &
        8.028_WP, 16.056_WP, 24.024_WP ]
    CHARACTER(len=6), PARAMETER :: shaft_keys(13) = [ &
        "lss_1p", "lss_2p", "lss_3p", "lss_6p", &
        "is1_1p","is1_2p","is1_3p", &
        "is2_1p","is2_2p","is2_3p", &
        "hss_1p", "hss_2p", "hss_3p" ]

    ! Gear mesh frequencies at min and max rpm
    REAL(WP), PARAMETER :: mesh_min(9) = [ &
        10.30_WP, 20.6_WP, 30.9_WP, &
        48.952_WP, 97.904_WP, 146.856_WP, &
        93.423_WP, 186.846_WP, 280.269_WP ]
    REAL(WP), PARAMETER :: mesh_max(9) = [ &
        16.48_WP, 32.96_WP, 49.44_WP, &
        78.3_WP, 156.6_WP, 234.9_WP, &
        149.492_WP, 298.984_WP, 448.476_WP ]
    CHARACTER(len=11), PARAMETER :: mesh_keys(9) = [ &
        "lss_mesh_1p","lss_mesh_2p","lss_mesh_3p", &
        "ims_mesh_1p","ims_mesh_2p","ims_mesh_3p", &
        "hss_mesh_1p","hss_mesh_2p","hss_mesh_3p" ]

    ! Base amplitudes (hardcoded)
    REAL(WP), PARAMETER :: A_base_shaft(13) = [ &
        0.6_WP, 0.3_WP, 1.0_WP, 0.2_WP, &
        0.4_WP, 0.2_WP, 0.6_WP, &
        0.3_WP, 0.15_WP, 0.4_WP, &
        0.2_WP, 0.1_WP, 0.3_WP ]
    REAL(WP), PARAMETER :: A_base_mesh(9) = [ &
        0.8_WP, 0.8_WP, 0.8_WP, &
        0.5_WP, 0.5_WP, 0.5_WP, &
        0.3_WP, 0.3_WP, 0.3_WP ]

    INTEGER(I32) :: n_total, i, m, pos_p, pos_underscore, n
    REAL(WP), ALLOCATABLE :: freqs(:), amps_source(:)
    REAL(WP) :: H, rms, fn, rpm_, damping_, alpha_mesh_, p_shaft_
    CHARACTER(len=20), ALLOCATABLE :: keys_work(:)
    CHARACTER(len=20) :: key_part

    
    ! Defaults
    rpm_        = 9.6_WP ; if (PRESENT(rpm))        rpm_        = rpm
    damping_    = 0.02_WP; if (PRESENT(damping))    damping_    = damping
    p_shaft_    = 1.0_WP ; if (PRESENT(p_shaft))    p_shaft_    = p_shaft
    alpha_mesh_ = 0.8_WP ; if (PRESENT(alpha_mesh)) alpha_mesh_ = alpha_mesh

    ! Clip rpm
    rpm_clipped = min(max(rpm_, min_rpm), max_rpm)
    factor = (rpm_clipped - min_rpm) / (max_rpm - min_rpm)

    ! Interpolate frequencies
    n_total = size(shaft_min) + size(mesh_min)
    ALLOCATE(freqs(n_total), amps_source(n_total), keys_work(n_total))

    ! Shaft frequencies and amplitudes
    do i = 1, size(shaft_min)
        freqs(i) = shaft_min(i) + factor * (shaft_max(i) - shaft_min(i))
        keys_work(i) = shaft_keys(i)

        pos_underscore = index(trim(shaft_keys(i)), "_", back=.true.)
        pos_p         = index(trim(shaft_keys(i)), "p", back=.true.)

        if (pos_underscore > 0 .and. pos_p > pos_underscore) then
            key_part = trim(shaft_keys(i)(pos_underscore+1:pos_p-1))
            read(key_part, *) n
        else
            n = 1
        end if

        amps_source(i) = A_base_shaft(i) * (1.0_WP / real(n, WP) ** p_shaft_)
    end do

    ! Mesh frequencies and amplitudes
    do i = 1, size(mesh_min)
        pos_underscore = index(trim(mesh_keys(i)), "_", back=.true.)
        pos_p         = index(trim(mesh_keys(i)), "p", back=.true.)

        if (pos_underscore > 0 .and. pos_p > pos_underscore) then
            key_part = trim(mesh_keys(i)(pos_underscore+1:pos_p-1))
            read(key_part, *) m
        else
            m = 1
        end if

        freqs(size(shaft_min) + i) = mesh_min(i) + factor * (mesh_max(i) - mesh_min(i))
        keys_work(size(shaft_min) + i) = mesh_keys(i)
        amps_source(size(shaft_min) + i) = A_base_mesh(i) * exp(-alpha_mesh_ * (real(m, WP) - 1.0_WP))
    end do

    ! Modal transfer function H(f)
    do i = 1, n_total
        modal_sum = 0.0_WP
        do m = 1, size(torsional_modes)
            fn = torsional_modes(m)
            modal_sum = modal_sum + 1.0_WP /((1.0_WP - (freqs(i)/fn)**2)**2 + (2.0_WP * damping_ * freqs(i)/fn)**2)
        end do
        H = sqrt(modal_sum)
        amps_source(i) = amps_source(i) * H
    end do

    ! Normalize to unit RMS: sqrt(sum(A^2/2)) = 1
    rms = sqrt(sum(amps_source**2)/2.0_WP)
    amps_source = amps_source / rms

    ! Build output array
    ALLOCATE(freqs_amp(n_total, 2))
    freqs_amp(:,1) = freqs
    freqs_amp(:,2) = amps_source
    ALLOCATE(keys(n_total))
    keys = keys_work

    DEALLOCATE(freqs, amps_source, keys_work)

END SUBROUTINE drivetrain10MW_excitation_spectrum


END MODULE TurbineTypes