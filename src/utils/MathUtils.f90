MODULE MathUtils

USE FFTW3
USE Kinds, ONLY: WP, I32, PI
USE, INTRINSIC :: iso_c_binding

IMPLICIT NONE
PUBLIC

INTERFACE
    SUBROUTINE qsort(base, nmemb, size, compar) BIND(c, name='qsort')
        IMPORT :: C_PTR, C_SIZE_T, C_FUNPTR
        TYPE(C_PTR), VALUE :: base
        INTEGER(C_SIZE_T), VALUE :: nmemb
        INTEGER(C_SIZE_T), VALUE :: size
        TYPE(C_FUNPTR), VALUE :: compar
    END SUBROUTINE qsort
END INTERFACE


CONTAINS


SUBROUTINE format_elapsed(start_time, elapsed_display, tag)
    ! Returns elapsed time since start_time and a human-friendly unit tag
    REAL(WP), INTENT(IN)  :: start_time
    REAL(WP), INTENT(OUT) :: elapsed_display
    CHARACTER(len=*), INTENT(OUT) :: tag

    REAL(WP) :: end_time, elapsed_time

    CALL cpu_time(end_time)
    elapsed_time = end_time - start_time

    elapsed_display = elapsed_time
    if (elapsed_display > 60.0_WP .and. elapsed_display <= 3600.0_WP) then
        elapsed_display = elapsed_display / 60.0_WP
        tag = "min"
    elseif (elapsed_display > 3600.0_WP .and. elapsed_display <= 86400.0_WP) then
        elapsed_display = elapsed_display / 3600.0_WP
        tag = "h"
    elseif (elapsed_display > 86400.0_WP) then
        elapsed_display = elapsed_display / 86400.0_WP
        tag = "days"
    else
        tag = "s"
    end if

END SUBROUTINE format_elapsed


FUNCTION remove_duplicate_nodes(x_all, delete_last) RESULT(mask)
    !> Removes duplicate nodes shared between consecutive structural members.
    !>
    !> Given an array of nodes of shape (Nm, Nn, 3), where Nm is the number of
    !> members and Nn the number of nodes per member, this function identifies
    !> nodes that are duplicated at the junction between member i and member i+1.
    !> The last node of member i is assumed to coincide with the first node of
    !> member i+1. One of them is marked for removal based on `delete_last`.
    !>
    !> Arguments:
    !>   x_all      : REAL(WP), dimension(:,:,:)  (Nm, Nn, 3)
    !>   delete_last: LOGICAL  If .true., the first node of the following member
    !>                          is removed; if .false., the last node of the
    !>                          previous member is removed.
    !>
    !> Returns:
    !>   mask       : LOGICAL, dimension(Nm*Nn)  Flat mask with .true. for nodes to keep
    !>     
    REAL(WP), INTENT(IN) :: x_all(:,:,:)   ! (Nm, Nn, 3)
    LOGICAL,  INTENT(IN) :: delete_last
    LOGICAL, ALLOCATABLE :: mask(:)

    integer(i32) :: Nm, Nn, i, idx_prev, idx_next
    REAL(WP) :: tol
    LOGICAL :: are_equal

    Nm = size(x_all, 1)
    Nn = size(x_all, 2)

    allocate(mask(Nm*Nn))
    mask = .true.   ! start keeping all nodes

    ! Tolerance for comparing coordinates (relative)
    tol = 1.0e-6_WP

    do i = 1, Nm-1
        ! Last node of current member
        idx_prev = (i-1)*Nn + Nn
        ! First node of next member
        idx_next = i*Nn + 1

        ! Compare coordinates (simple Euclidean distance, but we can check
        ! each coordinate difference)
        are_equal = all(abs(x_all(i, Nn, :) - x_all(i+1, 1, :)) < tol)

        if (are_equal) then
            if (delete_last) then
                ! Remove first node of next member
                mask(idx_next) = .false.
            else
                ! Remove last node of current member
                mask(idx_prev) = .false.
            end if
        end if
    end do

END FUNCTION remove_duplicate_nodes


FUNCTION divide_span(x) RESULT(dx)
    !> Computes representative span intervals (dx) for an ordered spatial vector.
    !> Equivalent to Python's divide_span.
    REAL(WP), INTENT(IN) :: x(:)
    REAL(WP), ALLOCATABLE :: dx(:)

    INTEGER(I32) :: n, i
    REAL(WP), ALLOCATABLE :: x_work(:), inf(:), sup(:)
    LOGICAL :: flipped

    n = size(x)
    allocate(x_work(n), inf(n), sup(n), dx(n))

    ! Check direction and possibly reverse
    flipped = .false.
    if (x(1) > x(n)) then
        x_work = x(n:1:-1)   ! reverse
        flipped = .true.
    else
        x_work = x
    end if

    ! Compute intervals
    do i = 1, n
        if (i == 1) then
            inf(i) = x_work(i)
        else
            inf(i) = 0.5_wp * (x_work(i) + x_work(i-1))
        end if

        if (i == n) then
            sup(i) = x_work(i)
        else
            sup(i) = 0.5_wp * (x_work(i) + x_work(i+1))
        end if

        dx(i) = sup(i) - inf(i)
    end do

    ! If we had reversed, reverse dx back
    if (flipped) then
        dx = dx(n:1:-1)
    end if

    ! Check consistency: all dx must have the same sign
    if (.not. (all(dx > 0.0_wp) .or. all(dx < 0.0_wp))) then
        error stop "divide_span: check input, non-monotonic grid"
    end if

END FUNCTION divide_span


SUBROUTINE compute_rfft(array, nt, dt, skipf, remove_zero, array_out, freqs)
    !> Computes the Real Fast Fourier Transform (RFFT) of a 3D time-series array.
    !> Equivalent to Python's compute_rfft.
    !>
    !> Arguments:
    !>   array  : real(wp), dimension(:,:,:)   Input array (nt, Nnodes, 3)
    !>   nt     : integer                      Number of time steps
    !>   dt     : real(wp)                     Time step size [s]
    !>   skipf  : integer, optional            Frequency decimation factor
    !>   remove_zero : logical, optional       Drop DC frequency (default .true.)
    !>
    !> Outputs:
    !>   f      : complex(wp), allocatable     Frequency-domain spectrum (nfreq, Nnodes, 3)
    !>   freqs  : real(wp), allocatable        Frequency bins [Hz]

    REAL(WP)                , INTENT(IN)           :: array(:,:,:)
    INTEGER(I32)            , INTENT(IN)           :: nt
    REAL(WP)                , INTENT(IN)           :: dt
    INTEGER(I32)            , INTENT(IN), OPTIONAL :: skipf
    LOGICAL                 , INTENT(IN), OPTIONAL :: remove_zero
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT)          :: array_out(:,:,:)
    REAL(WP), ALLOCATABLE   , INTENT(OUT)          :: freqs(:)

    ! Local variables
    INTEGER(I32) :: nfreq, n_signals, i, j, k, idx, pos
    INTEGER(I32) :: skipf_, nfreq_new, count_pos
    LOGICAL :: remove_zero_
    REAL(WP) :: win_norm
    REAL(WP), ALLOCATABLE :: window(:), mean_vals(:,:)
    REAL(WP), ALLOCATABLE :: data_in(:,:,:)
    REAL(WP), ALLOCATABLE :: data_flat(:,:)
    COMPLEX(WP), ALLOCATABLE :: out_flat(:,:)
    COMPLEX(WP), ALLOCATABLE :: f_temp(:,:,:)
    REAL(WP), ALLOCATABLE :: freqs_temp(:)
    COMPLEX(WP), ALLOCATABLE :: f_filt(:,:,:)
    REAL(WP), ALLOCATABLE :: freqs_filt(:)

    TYPE(C_PTR) :: plan
    INTEGER(C_INT) :: rank
    INTEGER(C_INT), ALLOCATABLE :: n_(:), inembed(:), onembed(:)
    INTEGER(C_INT) :: howmany, istride, idist, ostride, odist, flags


    ! Defaults
    skipf_ = 1           ; if (PRESENT(skipf))       skipf_       = skipf
    remove_zero_ = .true.; if (PRESENT(remove_zero)) remove_zero_ = remove_zero

    n_signals = size(array,2) * size(array,3)
    nfreq = nt/2 + 1

    ! Hanning window and its mean
    ALLOCATE(window(nt))
    do i = 1, nt
        window(i) = 0.5_WP - 0.5_WP*cos(2.0_WP*PI * REAL(i-1,WP) / REAL(nt-1, WP))
    end do
    win_norm = sum(window) / nt

    ! Substract mean along time axis
    ALLOCATE(mean_vals(size(array,2), size(array,3)))
    mean_vals = sum(array, dim=1)/nt
    ALLOCATE(data_in(nt, size(array,2), size(array,3)))
    do i = 1, nt
        data_in(i,:,:) = array(i,:,:) - mean_vals
    end do
    
    ! Apply window
    do i = 1, nt
        data_in(i,:,:) = data_in(i,:,:) * window(i)
    end do
    
    ! Flatten to (nt, n_signals) for FFTW many transforms
    ALLOCATE(data_flat(nt, n_signals))
    do k = 1, size(array,3)
        do j = 1, size(array,2)
            idx = (k-1)*size(array,2) + j
            data_flat(:,idx) = data_in(:,j,k)
        end do
    end do
    
    ! Prepare FFTW plan for many 1D real-to-complex transforms
    ALLOCATE(out_flat(nfreq, n_signals))
    rank = 1
    ALLOCATE(n_(1))
    n_(1) = nt
    ALLOCATE(inembed(1), onembed(1))
    inembed(1) = nt
    onembed(1) = nfreq
    howmany = n_signals
    istride = 1
    idist = nt
    ostride = 1
    odist = nfreq
    flags = FFTW_ESTIMATE

    plan = fftw_plan_many_dft_r2c(rank, n_, howmany, &
                                  data_flat, inembed, istride, idist, &
                                  out_flat, onembed, ostride, odist, flags)

    if (.not. c_associated(plan)) error stop "compute_rfft: FFTW plan creation failed"

    ! Execute FFT
    CALL fftw_execute_dft_r2c(plan, data_flat, out_flat)

    ! Destroy plan
    CALL fftw_destroy_plan(plan)

    ! Sacale and reshape output
    ALLOCATE(array_out(nfreq, size(array,2), size(array,3)))
    do k = 1, size(array,3)
        do j = 1, size(array,2)
            idx = (k-1) * size(array,2) + j
            array_out(:,j,k) = 2.0_WP * out_flat(:,idx) / (nt * win_norm)
        end do
    end do

    ! Frequencies
    ALLOCATE(freqs(nfreq))
    do i = 1, nfreq
        freqs(i) = REAL(I-1, WP) / (nt * dt)
    end do

    ! Applu skipf decimation
    if (skipf_ > 1) then
        print '(A,I6)', 'Number of frequencies before skipf: ', nfreq
        nfreq_new = (nfreq + skipf_ - 1) / skipf_
        ALLOCATE(f_temp(nfreq_new, size(array,2), size(array,3)))
        ALLOCATE(freqs_temp(nfreq_new))
        do i = 1, nfreq_new
            pos = (i-1) * skipf_ + 1
            freqs_temp(i) = freqs(pos)
            f_temp(i,:,:) = array_out(pos,:,:)
        end do
        DEALLOCATE(array_out, freqs)
        ALLOCATE(array_out(nfreq_new, size(array,2), size(array,3)))
        array_out = f_temp
        ALLOCATE(freqs(nfreq_new))
        freqs = freqs_temp
        DEALLOCATE(f_temp, freqs_temp)
    end if

    ! Remove zero frequency
    if (remove_zero_) then
        count_pos = count(freqs > 0.0_WP)
        ALLOCATE(f_filt(count_pos, size(array,2), size(array,3)))
        ALLOCATE(freqs_filt(count_pos))
        pos = 1
        do i = 1, size(freqs)
            if (freqs(i) > 0.0_wp) then
                freqs_filt(pos) = freqs(i)
                f_filt(pos,:,:) = array_out(i,:,:)
                pos = pos + 1
            end if
        end do
        DEALLOCATE(array_out, freqs)
        ALLOCATE(array_out(count_pos, size(array,2), size(array,3)))
        array_out = f_filt
        ALLOCATE(freqs(count_pos))
        freqs = freqs_filt
        DEALLOCATE(f_filt, freqs_filt)
    end if

    DEALLOCATE(window, mean_vals, data_in, data_flat, out_flat, n_, inembed, onembed)

END SUBROUTINE compute_rfft


SUBROUTINE generate_timeseries_banded_sines(peaks, keys, t , zeta, nfreq, seed, &
                                            used_freqs, fcut, a, freqs_out)
    !> Reconstructs a time-series signal from spectral peak distributions.
    !> Equivalent to Python's generate_timeseries_banded_sines.
    !>
    !> Arguments:
    !>   peaks     : real(wp), dimension(:,:)  (N,2) Col1: frequency [Hz], Col2: RMS amplitude
    !>   keys      : character(len=20), dimension(:)  Component labels
    !>   t         : real(wp), dimension(:)     Time array [s]
    !>   zeta      : real(wp)                   Modal damping ratio
    !>   nfreq     : integer                    Spectral lines per Gaussian band
    !>   seed      : integer                    Random seed (approximate)
    !>   used_freqs: logical                    If .true., output used frequencies
    !>   fcut      : real(wp)                   Cutoff frequency to force banded reconstruction
    !>
    !> Outputs:
    !>   a         : real(wp), allocatable      Synthesized time series
    !>   freqs_out : real(wp), allocatable      Sorted unique active frequencies (if used_freqs)
    USE, INTRINSIC :: iso_fortran_env

    REAL(WP)         , INTENT(IN) :: peaks(:,:)
    CHARACTER(len=20), INTENT(IN) :: keys(:)
    REAL(WP)         , INTENT(IN) :: t(:)
    REAL(WP)         , INTENT(IN), OPTIONAL :: zeta
    INTEGER(I32)     , INTENT(IN), OPTIONAL :: nfreq
    INTEGER(I32)     , INTENT(IN), OPTIONAL :: seed
    LOGICAL          , INTENT(IN), OPTIONAL :: used_freqs
    REAL(WP)         , INTENT(IN), OPTIONAL :: fcut
    REAL(WP), ALLOCATABLE, INTENT(OUT) :: a(:)
    REAL(WP), ALLOCATABLE, INTENT(OUT), OPTIONAL :: freqs_out(:)

    ! Local variables
    INTEGER(I32) :: i, j, seed_size, nt, n_peaks, n_active, nfreq_, seed_
    INTEGER(I32) :: max_freq_list, n_kept
    REAL(WP) :: f0, Arms, sigma, phase, rnd, zeta_, fcut_, sqrt2 = sqrt(2.0_WP)
    REAL(WP), ALLOCATABLE :: fk(:), g(:), Arms_k(:), phases(:)
    REAL(WP), ALLOCATABLE :: freq_list(:)
    INTEGER(I32), ALLOCATABLE :: seed_arr(:)
    LOGICAL :: used_freqs_
    CHARACTER(len=20) :: key_label

    ! Defaults
    zeta_       = 0.02_WP; if (PRESENT(zeta))       zeta_       = zeta
    nfreq_      = 50_I32 ; if (PRESENT(nfreq))      nfreq_      = nfreq
    seed_       = 42_I32 ; if (PRESENT(seed))       seed_       = seed
    fcut_       = 10.0_WP; if (PRESENT(fcut))       fcut_       = fcut
    used_freqs_ = .false.; if (PRESENT(used_freqs)) used_freqs_ = used_freqs
    
    nt = size(t)
    ALLOCATE(a(nt), source = 0.0_WP)
    n_peaks = size(peaks,1)
    if (size(keys) < n_peaks) error stop "generate_timeseries_banded_sines: keys shorter than peaks"
    if (size(peaks,2) < 2) error stop "generate_timeseries_banded_sines: peaks must have at least 2 columns"

    if (used_freqs_) then
        max_freq_list = max(1_I32, n_peaks * nfreq_)
        ALLOCATE(freq_list(max_freq_list))
        n_active = 0_I32
    else
        ALLOCATE(freq_list(0))
        n_active = 0_I32
    end if
    
    ! Initialize radnom seed
    CALL random_seed(size=seed_size)
    ALLOCATE(seed_arr(seed_size))
    CALL random_seed(put=seed_arr)
    
    do i = 1, n_peaks
        f0   = peaks(i,1)
        Arms = peaks(i,2)
        key_label = trim(keys(i))

        ! Gear mesh or high frequency -> banded Gaussian
        if (index(key_label, "mesh") > 0 .or. f0 > fcut_) then
            sigma = zeta_ * f0
            if (sigma <= 0.0_WP) cycle

            ALLOCATE(fk(nfreq_))
            do j = 1, nfreq_
                fk(j) = f0 - 4.0_WP * sigma + (j-1) * (8.0_WP * sigma) / REAL(nfreq_-1, WP)
            end do

            ALLOCATE(g(nfreq_))
            g = exp(-0.5_WP * ((fk-f0)/sigma)**2)
            g = g /sqrt(sum(g**2))

            ALLOCATE(Arms_k(nfreq_))
            Arms_k = Arms * g

            ALLOCATE(phases(nfreq_))
            do j = 1, nfreq_
                CALL random_number(rnd)
                phases(j) = 2.0_WP * PI * rnd
            end do

            ! Sum contribution
            do j = 1, nfreq_
                a = a + Arms_k(j) * sin(2.0_WP * PI * fk(j) * t + phases(j))
            end do

            if (used_freqs_) then
                do j = 1, nfreq_
                    n_active = n_active + 1
                    freq_list(n_active) = fk(j)
                end do
            end if

            DEALLOCATE (fk, g, Arms_k, phases)

        else

            ! Pure tone (shaft)
            CALL random_number(rnd)
            phase = 2.0_WP * PI * rnd
            a = a + sqrt2 * Arms * sin(2.0_WP * PI * f0 * t + phase)

            if (used_freqs_) then
                n_active = n_active + 1
                freq_list(n_active) = f0
            end if

        end if
    end do

    ! If used_freqs, sort and remove duplicates
    if (used_freqs_) then
        if (n_active > 0) then
            CALL sort_real_c(freq_list(1:n_active))
            ALLOCATE(freqs_out(n_active))
            n_kept = 0_I32
            do i = 1, n_active
                if (i == 1) then
                    n_kept            = n_kept + 1
                    freqs_out(n_kept) = freq_list(i)
                else if ( abs(freq_list(i) - freq_list(i-1)) > 1.0e-6_WP ) then
                    n_kept            = n_kept + 1
                    freqs_out(n_kept) = freq_list(i)
                end if
            end do
            if (n_kept < n_active) freqs_out = freqs_out(1:n_kept)
        else
             ALLOCATE(freqs_out(0))
        end if
    else
        if (PRESENT(freqs_out)) ALLOCATE (freqs_out(0))
    end  if

    DEALLOCATE(freq_list, seed_arr)

END SUBROUTINE generate_timeseries_banded_sines


FUNCTION filter_non_usefull_freqs(freqs, mantain_freqs, freqs_over, freqs_under) RESULT(mask)
    !> Generates a frequency-mask to filter out non-essential spectral bins.
    !> Equivalent to Python's filter_non_usefull_freqs.
    REAL(WP), INTENT(IN) :: freqs(:)
    REAL(WP), INTENT(IN) :: mantain_freqs(:)
    REAL(WP), INTENT(IN), OPTIONAL :: freqs_over
    REAL(WP), INTENT(IN), OPTIONAL :: freqs_under
    LOGICAL, ALLOCATABLE :: mask(:)

    ! Local variables
    INTEGER(I32) :: n, i, j
    REAL(WP) :: df, freqs_over_, freqs_under_
    LOGICAL :: equally_spaced
    REAL(WP), ALLOCATABLE :: diffs(:)
    LOGICAL, ALLOCATABLE :: filter_zone(:)


    ! Defaults
    freqs_over_  = minval(freqs); if (PRESENT(freqs_over))  freqs_over_  = freqs_over
    freqs_under_ = maxval(freqs); if (PRESENT(freqs_under)) freqs_under_ = freqs_under

    n = size(freqs)
    ALLOCATE(mask(n))
    mask = .true.

    ! Compute df assuming unifrom spacing
    if (n > 1) then
        df = freqs(2) - freqs(1)
        ! Check equal spacing
        ALLOCATE(diffs(n-1))
        diffs = freqs(2:n) - freqs(1:n-1)
        equally_spaced = all(abs(diffs-df) < 1.0e-6_WP)
        DEALLOCATE(diffs)
        if (.not. equally_spaced) error stop "filter_non_usefull_freqs: freqs must be equally spaced."
    else
        df = 0.0_WP
    end if

    ! Build filter zone
    ALLOCATE(filter_zone(n))
    filter_zone = .true.
    if (PRESENT(freqs_over))  filter_zone = filter_zone .and. (freqs >= freqs_over_)
    if (PRESENT(freqs_under)) filter_zone = filter_zone .and. (freqs <= freqs_under_)

    ! In filter zone, initially set mask to false
    where(filter_zone) mask = .false.

    ! For each mantained frequency, set mask true if within df
    do i = 1, n
        if (.not. filter_zone(i)) cycle
        do j = 1, size(mantain_freqs)
            if (abs(freqs(i) - mantain_freqs(j)) <= df) then
                mask(i) = .true.
                exit
            end if
        end do
    end do

    DEALLOCATE(filter_zone)

END FUNCTION filter_non_usefull_freqs


FUNCTION compare_real(a,b) BIND(c) RESULT(res)
    !> Comparison function for qsort: ascending order for real(c_double)

    TYPE(C_PTR), VALUE :: a, b
    INTEGER(C_INT) :: res
    REAL(C_DOUBLE), POINTER :: fa, fb

    CALL c_f_pointer(a, fa)
    CALL c_F_pointer(b, fb)

    if (fa < fb) then
        res = -1_C_INT
    else if (fa > fb) then
        res = 1_C_INT
    else
        res = 0_C_INT
    end if

END FUNCTION compare_real


SUBROUTINE sort_real_c(arr)
    !> Sorts a real(wp) array using C's qsort (ascending).

    REAL(WP), INTENT(INOUT), TARGET :: arr(:)
    REAL(WP), POINTER :: parr(:)
    INTEGER(C_SIZE_T) :: nmemb, size_elem
    TYPE(C_PTR) :: ptr

    parr => arr
    nmemb = size(parr, kind=C_SIZE_T)
    size_elem = c_sizeof(parr(1))
    ptr = c_loc(parr(1))

    CALL qsort(ptr, nmemb, size_elem, c_funloc(compare_real))
END SUBROUTINE sort_real_c


END MODULE MathUtils