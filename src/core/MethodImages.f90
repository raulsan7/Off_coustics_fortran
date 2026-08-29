MODULE MethodImages

USE omp_lib
USE Kinds, ONLY: WP, I32
USE WindTurbine, ONLY: WindTurbine_t
USE AcousticSolver, ONLY: AcousticSolver_t

IMPLICIT NONE

PRIVATE
PUBLIC :: MethodImages_t, WindTurbine_t


! ------------------------
! Derived Types Definition
! ------------------------
TYPE :: ImageSystem_t

    REAL(WP)    , ALLOCATABLE :: nodes_pos(:,:)     ! [m] Nodes coordinates (Nnodes, 3)
    COMPLEX(WP) , ALLOCATABLE :: force(:,:,:)       ! [N] Force array (Nfreqs, Nnodes, 3)
    INTEGER(I32), ALLOCATABLE :: BC_all(:)          ! [-] Boundary conditions array

END TYPE ImageSystem_t


TYPE, EXTENDS(AcousticSolver_t) :: MethodImages_t
    
    INTEGER(I32) :: N_images         = 30           ! [-] Number of image levels to compute
    INTEGER(I32) :: default_N_images = 30
    INTEGER(I32) :: Upper_BC         = -1           ! [-] Upper Boundary condition -1 for  p = 0
    INTEGER(I32) :: Lower_BC         = +1           ! [-] Lower Boundary condition +1 for dp = 0
    REAL(WP)     :: Upper_HBC        = 0.0_WP       ! [m] z-coordinate of the upper reflecting plane (e.g. free surface)
    REAL(WP)     :: Lower_HBC        = -30.0        ! [m] z-coordinate of the lower reflecting plane (e.g. seabed)
    REAL(WP)     :: up_R             = 1.0_WP       ! [-] Attenuation coefficient for upper boundary reflections
    REAL(WP)     :: lw_R             = 0.5_WP       ! [-] Attenuation coefficient for lower boundary reflections

    CLASS(WindTurbine_t), ALLOCATABLE :: turbines(:) 
    TYPE(ImageSystem_t) , ALLOCATABLE :: image_system(:)

    CONTAINS
    ! --- Public type-bound procedures --- !
    GENERIC   :: init        => init_array, init_single
    PROCEDURE :: init_array  => init_MethodImages
    PROCEDURE :: init_single => init_MethodImages_single
    PROCEDURE :: set_N_images
    PROCEDURE :: build_all_image_system
    PROCEDURE :: build_image_system
    PROCEDURE :: restore_default_images
    PROCEDURE :: method_of_images
    PROCEDURE :: dipole_pressure_images

    ! --- Deferred(abstract) procedures --- !
    PROCEDURE :: get_name           => get_name_MethodImages
    PROCEDURE :: compute_pressure   => compute_pressure_MethodImages

END TYPE MethodImages_t


CONTAINS

! ---------- DEFERRED ---------- !
SUBROUTINE init_MethodImages(self, turbines, N_images, c_wat, rho_wat, &
                             Upper_BC, Lower_BC, Upper_HBC, Lower_HBC, &
                             attenuation_upper, attenuation_lower, eps, p_ref, cluster)

    CLASS(MethodImages_t), INTENT(INOUT)        :: self
    CLASS(WindTurbine_t) , INTENT(IN), OPTIONAL :: turbines(:)
    INTEGER(I32), INTENT(IN), OPTIONAL          :: N_images
    INTEGER(I32), INTENT(IN), OPTIONAL          :: Upper_BC, Lower_BC
    REAL(WP)    , INTENT(IN), OPTIONAL          :: c_wat, rho_wat
    REAL(WP)    , INTENT(IN), OPTIONAL          :: Upper_HBC, Lower_HBC
    REAL(WP)    , INTENT(IN), OPTIONAL          :: attenuation_upper, attenuation_lower
    REAL(WP)    , INTENT(IN), OPTIONAL          :: eps, p_ref
    LOGICAL     , INTENT(IN), OPTIONAL          :: cluster


    ! Defaults
    if (PRESENT(N_images))          self % N_images  = N_images
    if (PRESENT(c_wat))             self % c         = c_wat
    if (PRESENT(rho_wat))           self % rho       = rho_wat
    if (PRESENT(Upper_BC))          self % Upper_BC  = Upper_BC
    if (PRESENT(Lower_BC))          self % Lower_BC  = Lower_BC
    if (PRESENT(Upper_HBC))         self % Upper_HBC = Upper_HBC
    if (PRESENT(Lower_HBC))         self % Lower_HBC = Lower_HBC
    if (PRESENT(attenuation_upper)) self % up_R      = attenuation_upper
    if (PRESENT(attenuation_lower)) self % lw_R      = attenuation_lower
    if (PRESENT(eps))               self % eps       = eps
    if (PRESENT(p_ref))             self % p_ref     = p_ref
    if (PRESENT(cluster))           self % cluster   = cluster

    
    self % default_N_images = self % N_images
    
    if (self % Upper_HBC <= self % Lower_HBC) error stop "MethodImages % init: Upper_HBC must be strictly greater than Lower_HBC"
    
    if (PRESENT(turbines)) then
        if (ALLOCATED(self % turbines)) DEALLOCATE(self % turbines)
        ALLOCATE(self % turbines, source = turbines)
    end if

    CALL self % build_all_image_system()

END SUBROUTINE init_MethodImages


SUBROUTINE init_MethodImages_single(self, turbine, N_images, c_wat, rho_wat, &
                                    Upper_BC, Lower_BC, Upper_HBC, Lower_HBC, &
                                    attenuation_upper, attenuation_lower, eps, p_ref, cluster)

    CLASS(MethodImages_t), INTENT(INOUT) :: self
    CLASS(WindTurbine_t) , INTENT(IN)    :: turbine
    INTEGER(I32), INTENT(IN), OPTIONAL   :: N_images
    REAL(WP)    , INTENT(IN), OPTIONAL   :: c_wat, rho_wat
    REAL(WP)    , INTENT(IN), OPTIONAL   :: Upper_HBC, Lower_HBC
    INTEGER(I32), INTENT(IN), OPTIONAL   :: Upper_BC, Lower_BC
    REAL(WP)    , INTENT(IN), OPTIONAL   :: attenuation_upper, attenuation_lower
    REAL(WP)    , INTENT(IN), OPTIONAL   :: eps, p_ref
    LOGICAL     , INTENT(IN), OPTIONAL   :: cluster


    if (allocated(self%turbines)) deallocate(self%turbines)
    allocate(self%turbines(1), source=turbine)

    ! Delegate the rest of initialization to the array-based init
    call init_MethodImages(self, N_images=N_images, c_wat=c_wat, rho_wat=rho_wat, &
                           Upper_BC=Upper_BC, Lower_BC=Lower_BC, Upper_HBC=Upper_HBC, Lower_HBC=Lower_HBC, &
                           attenuation_upper=attenuation_upper, attenuation_lower=attenuation_lower, &
                           eps=eps, p_ref=p_ref, cluster=cluster)

END SUBROUTINE init_MethodImages_single


FUNCTION get_name_MethodImages(self) RESULT(name)
    CLASS(MethodImages_t), INTENT(INOUT) :: self
    CHARACTER(len=:)     , ALLOCATABLE   :: name
    name = "Images Method"; self % solver_name = trim(name)
END FUNCTION get_name_MethodImages


SUBROUTINE compute_pressure_MethodImages(self, observers, block_size, total_pressure)
        CLASS(MethodImages_t)   , INTENT(INOUT)        :: self
        REAL(WP)                , INTENT(IN)           :: observers(:,:)        ! [m] Observers coords (Nobs,3)
        INTEGER(I32)            , INTENT(IN), OPTIONAL :: block_size            ! [Pa] Pressure field (Nfreqs, Nobs)
        COMPLEX(WP), ALLOCATABLE, INTENT(OUT)          :: total_pressure(:,:)   ! [-] Number of observers to process simultaneously

    print*, self % cluster, observers, block_size, total_pressure

END SUBROUTINE compute_pressure_MethodImages


! ---------- HELPERS ---------- !
SUBROUTINE build_all_image_system(self)

    CLASS(MethodImages_t), INTENT(INOUT) :: self
    INTEGER(I32) :: i, num_turbines
    COMPLEX(WP), ALLOCATABLE :: F_corr(:,:,:)

    if (.not. ALLOCATED(self % turbines)) return
    num_turbines = size(self % turbines)
    
    if (ALLOCATED(self % image_system)) deallocate(self % image_system)
    allocate(self % image_system(num_turbines))
    
    do i = 1, num_turbines
        F_corr = self % turbines(i) % get_impedance_corrected_force(self % c)
        CALL self % build_image_system(self % turbines(i) % x, F_corr, self % image_system(i))
    end do

END SUBROUTINE build_all_image_system


SUBROUTINE build_image_system(self, nodes_pos_real, force, img_sys)
    CLASS(MethodImages_t), INTENT(IN) :: self
    REAL(WP), INTENT(IN)              :: nodes_pos_real(:,:)   ! (Nnodes, 3)
    COMPLEX(WP), INTENT(IN)           :: force(:,:,:)          ! (Nfreqs, Nnodes, 3) or (Nfreqs, 3, Nnodes)
    TYPE(ImageSystem_t), INTENT(OUT)  :: img_sys

    REAL(WP), ALLOCATABLE :: zi(:), z_all(:)
    COMPLEX(WP), ALLOCATABLE :: force_prep(:,:,:), Force_out(:,:,:)
    INTEGER(I32), ALLOCATABLE :: parent(:), BC_all(:)
    INTEGER(I32) :: Nnodes, Nfreqs, total_nodes, j

    Nnodes = size(nodes_pos_real, 1)

    ! Ensure shape standard: (Nfreqs, 3, Nnodes)
    if (size(force, 3) == 3 .and. size(force, 2) == Nnodes) then
        Nfreqs = size(force, 1)
        allocate(force_prep(Nfreqs, 3, Nnodes))
        do j = 1, Nnodes
            force_prep(:, 1, j) = force(:, j, 1)
            force_prep(:, 2, j) = force(:, j, 2)
            force_prep(:, 3, j) = force(:, j, 3)
        end do
    else
        force_prep = force
    end if

    zi = nodes_pos_real(:, 3)

    CALL self % method_of_images(zi, force_prep, z_all, Force_out, parent, BC_all)

    total_nodes = size(z_all)
    allocate(img_sys%nodes_pos(3, total_nodes))
    allocate(img_sys%force(size(Force_out, 1), 3, total_nodes))
    allocate(img_sys%BC_all(total_nodes))

    do j = 1, total_nodes
        img_sys%nodes_pos(1, j) = nodes_pos_real(parent(j), 1)
        img_sys%nodes_pos(2, j) = nodes_pos_real(parent(j), 2)
        img_sys%nodes_pos(3, j) = z_all(j)
    end do

    img_sys%force  = Force_out
    img_sys%BC_all = BC_all

END SUBROUTINE build_image_system


SUBROUTINE set_N_images(self, new_N_images)
    CLASS(MethodImages_t), INTENT(INOUT) :: self
    INTEGER(I32), INTENT(IN)             :: new_N_images

    self % N_images = new_N_images
    CALL self % build_all_image_system()
END SUBROUTINE set_N_images


SUBROUTINE restore_default_images(self)
    CLASS(MethodImages_t), INTENT(INOUT) :: self

    if (self%N_images /= self%default_N_images) then
        self%N_images = self%default_N_images
        CALL self%build_all_image_system()
    end if
END SUBROUTINE restore_default_images


! ---------- METHOD ---------- !
SUBROUTINE method_of_images(self, zi, Force, z_all, Force_out, parent, BC_all)
    CLASS(MethodImages_t), INTENT(IN)         :: self
    REAL(WP), INTENT(IN)                      :: zi(:)
    COMPLEX(WP), INTENT(IN)                   :: Force(:,:,:)
    REAL(WP), ALLOCATABLE, INTENT(OUT)        :: z_all(:)
    COMPLEX(WP), ALLOCATABLE, INTENT(OUT)     :: Force_out(:,:,:)
    INTEGER(I32), ALLOCATABLE, INTENT(OUT)    :: parent(:)
    INTEGER(I32), ALLOCATABLE, INTENT(OUT)    :: BC_all(:)

    ! Local variables
    INTEGER(I32) :: Nnodes, Nfreqs, total_nodes, idx, i_node, img_lvl
    REAL(WP) :: z_upper, z_lower
    COMPLEX(WP), ALLOCATABLE :: F_upper(:,:), F_lower(:,:)
    CHARACTER(len=5) :: last_plane_upper, last_plane_lower
    INTEGER(I32) :: BC_u, BC_l


    Nnodes = size(zi); Nfreqs = size(Force, 1)
    total_nodes = Nnodes * (1 + 2 * self % N_images)

    ALLOCATE(z_all(total_nodes))
    ALLOCATE(Force_out(Nfreqs, 3, total_nodes))
    ALLOCATE(parent(total_nodes))
    ALLOCATE(BC_all(total_nodes))
    ALLOCATE(F_upper(Nfreqs, 3))
    ALLOCATE(F_lower(Nfreqs, 3))

    idx = 1
    do i_node = 1, Nnodes
        ! Real dipole nodes
        z_all(idx)         = zi(i_node)
        Force_out(:,:,idx) = Force(:,:,i_node)
        BC_all(idx)        = 1
        parent(idx)        = i_node
        idx                = idx + 1

        z_upper = zi(i_node)
        z_lower = zi(i_node)
        F_upper = Force(:,:,i_node)
        F_lower = Force(:,:,i_node)

        last_plane_upper = "upper"
        last_plane_lower = "lower"

        do img_lvl = 1, self % N_images
            ! ---------- Upper Reflection Chain ---------- !
            if (last_plane_upper == "upper") then
                z_upper = 2.0_WP * self % Lower_HBC - z_upper
                BC_u = self % Lower_BC
                last_plane_upper = "lower"
            else
                z_upper = 2.0_WP * self % Upper_HBC - z_upper
                BC_u = self % Upper_BC
                last_plane_upper = "upper"
            end if

            F_upper(:,3)       = - 1.0_WP * F_upper(:,3)
            F_upper            = F_upper * self % up_R

            z_all(idx)         = z_upper
            Force_out(:,:,idx) = F_upper
            BC_all(idx)        = BC_u
            parent(idx)        = i_node
            idx                = idx + 1

            ! ---------- Lower Reflection Chain ---------- !
            if (last_plane_lower == "lower") then
                z_lower = 2.0_WP * self % Upper_HBC - z_lower
                BC_l    = self % Upper_BC
                last_plane_lower = "upper"
            else 
                z_lower = 2.0_WP * self % Lower_HBC - z_lower
                BC_l    = self % Lower_BC
                last_plane_lower = "lower"
            end if

            F_lower(:,3) = -1.0_WP * F_lower(:,3)
            F_lower      = F_lower * self % lw_R

            z_all(idx)         = z_lower
            Force_out(:,:,idx) = F_lower
            BC_all(idx)        = BC_l
            parent(idx)        = i_node
            idx                = idx + 1
        end do
    end do


END SUBROUTINE method_of_images


SUBROUTINE dipole_pressure_images(self,obs_pos, freqs, nodes_pos, force, BC_all, p_out)
    USE Kinds, ONLY: I1, PI
    CLASS(MethodImages_t), INTENT(IN)         :: self
     REAL(WP), INTENT(IN)                  :: obs_pos(3)
     REAL(WP), INTENT(IN)                  :: freqs(:)
     REAL(WP), INTENT(IN)                  :: nodes_pos(:,:)     ! (Ntotal, 3)
     COMPLEX(WP), INTENT(IN)               :: force(:,:,:)       ! (Nfreqs, Ntotal, 3)
     INTEGER(I32), INTENT(IN)              :: BC_all(:)          ! (Ntotal)
     COMPLEX(WP), ALLOCATABLE, INTENT(OUT) :: p_out(:)           ! (Nfreqs)

    ! Local variables
    INTEGER(I32) :: Nfreqs, Ntotal, f, j
    REAL(WP) :: rx, ry, rz, r, r_inv, k
    COMPLEX(WP) :: F_dot_r, term1, green


    Nfreqs = size(freqs); Ntotal = size(nodes_pos, 1)
    ALLOCATE(p_out(Nfreqs))
    p_out = (0.0_WP, 0.0_WP)

    do j = 1, Ntotal
        rx = obs_pos(1) - nodes_pos(j,1)
        rx = obs_pos(2) - nodes_pos(j,2)
        rx = obs_pos(3) - nodes_pos(j,3)

        r = sqrt(rx*rx + ry*ry + rz*rz)
        if (r < self % eps) r = self % eps
        r_inv = 1.0_WP / r

        do f = 1, Nfreqs
            k = 2.0_WP * PI * freqs(f) / self % c

            F_dot_r = (force(f,j,1) * rx + force(f,j,2) + force(f,j,3)*rz) * r_inv
            term1   = -I1 * k + r_inv
            green   = exp(CMPLX(0.0_WP, k * r, kind=WP)) / (4.0_WP * PI *r)

            p_out(f) = p_out(f) + F_dot_r * term1 * green * REAL(BC_all(j), WP)
        end do
    end do




END SUBROUTINE dipole_pressure_images

END MODULE MethodImages