!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Green's functions.        !
!! Copyright (C) 2012                                                       !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Function calculation                             !
!!                                                                          !
!! Developers: Alessandro Pecchia, Gabriele Penazzi                         !
!! Former Conctributors: Luca Latessa, Aldo Di Carlo                        !
!!                                                                          !
!! libNEGF is free software: you can redistribute it and/or modify          !
!! it under the terms of the GNU Lesse General Public License as published  !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !
!!--------------------------------------------------------------------------!


module lib_param

  use ln_precision, only : dp
  use globals
  use mat_def
  use ln_structure, only : TStruct_info, TBasisCenters, TNeighbourMap, print_Tstruct
  use input_output
  use elph, only : Telph
  use phph
  use energy_mesh, only : mesh
  use interactions, only : TInteraction, TInteractionList
  use elphdd, only : ElPhonDephD, ElPhonDephD_create, ElPhonDephD_init
  use elphdb, only : ElPhonDephB, ElPhonDephB_create, ElPhonDephB_init
  use elphds, only : ElPhonDephS, ElPhonDephS_create, ElPhonDephS_init
  use elphinel, only : ElPhonInel, ElPhonInel_create, ElPhonInel_init
  use scba
  use ln_cache
#:if defined("MPI")
  use libmpifx_module, only : mpifx_comm
#:endif
  implicit none
  private

  public :: Tnegf, intArray, TEnGrid, TInteractionList

  public :: set_bp_dephasing
  public :: set_elph_dephasing
  public :: set_elph_block_dephasing
  public :: set_elph_s_dephasing
  public :: set_elph_inelastic
  public :: set_phph
  integer, public, parameter :: MAXNCONT=10

  type intArray
    integer, dimension(:), allocatable :: indexes
  end type intArray


  !! Structure used to define energy points for the integration
  !! For every point we define
  !!     path (1,2 or 3): the energy point belongs to a real axis
  !!     integration (1), a complex plane integration (2) or a
  !!     pole summation (3)
  !!     pt_path: relative point number within a single path
  !!     pt: absolute point number along the whole integration path
  !!     cpu: cpu assigned to the calculation of the given energy point
  !!     Ec: energy value
  !!     wght: a weight used in final summation to evaluate integrals
  type TEnGrid
     integer :: path
     integer :: pt_path
     integer :: pt
     integer :: cpu
     complex(dp) :: Ec
     complex(dp) :: wght
  end type TEnGrid

  !-----------------------------------------------------------------------------
  ! Buttiker probe dephasing
  type :: Tdeph_bp
     real(dp), allocatable, dimension(:) :: coupling
  end type Tdeph_bp

  !-----------------------------------------------------------------------------
  ! Structure to hold contact and reservoir data
  type Tcontact
    character(132) :: name
    logical :: tWriteSelfEnergy = .false.
    logical :: tReadSelfEnergy = .false.
    logical :: tWriteSurfaceGF = .false.
    logical :: tReadSurfaceGF = .false.
    complex(dp), dimension(:,:,:), allocatable :: SelfEnergy ! Electrode Self-Energy
    complex(dp), dimension(:,:,:), allocatable :: SurfaceGF  ! Electrode Surface Green Function
    real(dp) :: mu_n        ! Electrochemical potential (el)
    real(dp) :: mu_p        ! Electrochemical potential (hl)
    real(dp) :: mu          ! Electrochemical Potential (dft calculation)
    real(dp) :: contact_DOS ! Ficticious contact DOS
    logical  :: FictCont    ! Ficticious contact
    real(dp) :: kbT_dm      ! Electronic temperature
    real(dp) :: kbT_t       ! Electronic temperature
    type(z_DNS) :: HC       ! Contact Hamiltonian
    type(z_DNS) :: SC       ! Contact Overlap
    type(z_DNS) :: HMC      ! Device-Contact Hamiltonian
    type(z_DNS) :: SMC      ! Device-Contact Overlap
  end type Tcontact

  !> General libnegf container
  !! Contains input data, runtime quantities and output data
  type Tnegf
    !! Input parameters: set by library user
    !! General
    integer :: verbose
#:if defined("MPI")
    type(mpifx_comm) :: globalComm
    type(mpifx_comm) :: cartComm
    type(mpifx_comm) :: energyComm
    type(mpifx_comm) :: kComm
#:endif
    integer  :: ReadoldDM_SGFs    ! 0: Read 1: compute 2: comp & save
    integer  :: ReadoldT_SGFs     ! 0: Read 1: compute 2: comp & save
    character(len=LST) :: scratch_path    ! Folder for scratch work
    character(len=LST) :: out_path        ! Folder for output data
    real(dp) :: g_spin            ! spin degeneracy
    real(dp) :: delta             ! delta for G.F.
    real(dp) :: dos_delta         ! additional delta to force more broadening in the DOS
    integer  :: deltaModel        ! Used for phonon G.F. ! delta**2, 2*delta*w, Mingo's
    real(dp) :: wmax              ! Maximum frequency in Mingo's model
                                  ! See 'Numerical Heat Transfer, Part B', 51:333, 2007
    real(dp) :: eneconv           ! Energy conversion factor
    integer  :: spin = 1          ! spin component
    character(1) :: DorE          ! Density or En.Density

    !! Contacts info
    type(Tcontact), dimension(:), allocatable :: cont

    !! Contour integral
    integer :: Np_n(2)            ! Number of points for n
    integer :: Np_p(2)            ! Number of points for p
    integer :: Np_real            ! Number of points for integration over real axis
    integer :: n_kt               ! Number of kT extending integrations
    integer :: n_poles            ! Number of poles
    real(dp) :: Ec                ! conduction band edge
    real(dp) :: Ev                ! valence band edge

    !! Real axis
    real(dp) :: Emin              ! Tunneling or dos interval
    real(dp) :: Emax              !
    real(dp) :: Estep             ! Tunneling or dos E step

    !! Emitter and collector for transmission or Meir-Wingreen
    !! (only emitter in this case)
    integer, allocatable :: ni(:) ! ni: emitter contact list
    integer, allocatable :: nf(:) ! nf: collector contact list

    integer :: ndos_proj              ! number of LDOS interval ranges
    type(intArray), dimension(:), allocatable :: DOS_proj ! PDOS descriptor
                                                          ! contain matrix index
                                                          ! for PDOS projection
    type(intArray) :: TUN_proj    ! Array of TUN projection indices

    real(dp) :: DeltaEc           ! safe guard energy below Ec
    real(dp) :: DeltaEv           ! safe guard energy above Ev

    !! Runtime variables: used internally by the library
    type(format) :: form          ! Form of file-Hamiltonian
    logical  :: dumpHS            ! Used for debug
    real(dp) :: muref             ! reference elec.chem potential
    real(dp) :: E                 ! Holding variable
    real(dp) :: dos               ! Holding variable
    integer :: activecont         ! contact selfenergy
    integer :: min_or_max         ! in input: 0 take minimum, 1 take maximum mu
    integer :: refcont            ! reference contact (for non equilib)
    integer :: outer              ! flag switching computation of
                                  ! the Device/Contact DM
                                  ! 0 none; 1 upper block; 2 all

    real(dp) :: mu                !    chem potential used without contacts
    real(dp) :: mu_n              ! el chem potential used without contacts
    real(dp) :: mu_p              ! hl chem potential used without contacts
    real(dp) :: kbT               ! temperature used without contacts

    !! Note: H,S are partitioned immediately after input, therefore they are
    !! built runtime from input variable
    type(z_CSR), pointer :: H => null()    ! Points to externally allocated H
    type(z_CSR), pointer :: S => null()
    type(z_CSR), pointer :: rho => null()      ! Holding output Matrix
    type(z_CSR), pointer :: rho_eps => null()  ! Holding output Matrix
    logical    :: isSid           ! True if overlap S == Id
    logical    :: intHS           ! tells HS are internally allocated
    logical    :: intDM           ! tells DM is internally allocated

    type(TStruct_Info) :: str     ! system structure
    type(TBasisCenters) :: basis  ! local basis centers
    type(TNeighbourMap), dimension(:), allocatable :: neighbour_map
    integer :: iE                 ! Currently processed En point (index)
    complex(dp) :: Epnt           ! Currently processed En point (value)
    real(dp) :: kwght             ! currently processed k-point weight
    integer :: iKpoint            ! currently processed k-point (index)

    ! Energy grid information
    type(TEnGrid), dimension(:), allocatable :: en_grid
    integer :: local_en_points    ! Local number of energy points

    ! Array to store all k-points and k-weights
    real(dp), allocatable, dimension(:,:) :: kpoints
    real(dp), allocatable, dimension(:) :: kweights
    ! Array of local k-point indices
    integer, allocatable, dimension(:) :: local_k_index

    type(mesh) :: emesh           ! energy mesh for adaptive Simpson
    real(dp) :: int_acc           ! adaptive integration accuracy

    type(Telph) :: elph           ! electron-phonon data
    type(Tphph) :: phph           ! phonon-phonon data

    ! Many Body Interactions as array of pointers
    type(TInteractionList)  :: interactList
    type(TScbaDriver) :: scbaDriver

    !! Output variables: these arrays are filled by internal subroutines to store
    !! library outputs
    real(dp), dimension(:,:), allocatable :: tunn_mat
    real(dp), dimension(:,:), allocatable :: curr_mat
    real(dp), dimension(:,:), allocatable :: ldos_mat
    real(dp), dimension(:), allocatable :: currents

    ! These variables need to be done to clean up
    logical :: tOrthonormal = .false.
    logical :: tOrthonormalDevice = .false.
    integer :: numStates = 0
    logical :: tDestroyGr = .true.
    logical :: tDestroyGn = .true.
    logical :: tDestroyESH = .true.
    ! Buttiker Probes dephasing
    logical :: tDephasingBP = .false.
    logical :: tDephasingVE = .false.
    type(Tdeph_bp) :: bp_deph

    ! internal use only
    integer :: readOldSGF

    ! Work variable: surface green cache.
    class(TMatrixCache), allocatable :: surface_green_cache
    class(TMatrixCache), allocatable :: ESH
    ! These are pointers so they can be passed to inelastic
    class(TMatrixCache), pointer :: G_r => null()
    class(TMatrixCache), pointer :: G_n => null()

    contains

    procedure :: set_defaults => set_defaults
    procedure :: print_all_vars => print_all_vars
    procedure :: destroy_interactions => destroy_interactions

  end type Tnegf

contains

  !> Set buttiker probe dephasing
  subroutine set_bp_dephasing(negf, coupling)
    type(Tnegf), intent(inout) :: negf
    real(dp),  dimension(:), intent(in) :: coupling

    if (.not.allocated(negf%bp_deph%coupling)) then
       allocate(negf%bp_deph%coupling(size(coupling)))
    end if
    negf%bp_deph%coupling = coupling

    ! Initialize the cache space
    call initialize_cache_space(negf,'G_r')
    call initialize_cache_space(negf,'G_n')

  end subroutine set_bp_dephasing

  !> Set values for the local electron phonon dephasing model
  !! (elastic scattering only)
  subroutine set_elph_dephasing(negf, coupling, niter)
    type(Tnegf), intent(inout) :: negf
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    integer, intent(in) :: niter

    type(TInteractionNode), pointer :: node
    
    call negf%interactList%add(node)
    call elphondephd_create(node%inter)
    select type(pInter => node%inter)
    type is(ElPhonDephD)
      call elphondephd_init(pInter, negf%str, coupling, niter)
    class default
      stop 'ERROR: error of type downcast to ElPhonDephD'
    end select

    ! Initialize the cache space
    call initialize_cache_space(negf,'G_r')
    call initialize_cache_space(negf,'G_n')

  end subroutine set_elph_dephasing

  !> Set values for the semi-local electron phonon dephasing model
  !! (elastic scattering only)
  subroutine set_elph_block_dephasing(negf, coupling, orbsperatom, niter)
    type(Tnegf), intent(inout) :: negf
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    integer,  dimension(:), allocatable, intent(in) :: orbsperatom
    integer, intent(in) :: niter

    type(TInteractionNode), pointer :: node
    
    call negf%interactList%add(node)
    call elphondephd_create(node%inter)
    select type(pInter => node%inter)
    type is(ElPhonDephB)
      call elphondephb_init(pInter, negf%str, coupling, orbsperatom, niter)
    class default
      stop 'ERROR: error of type downcast to ElPhonDephB'
    end select

    ! Initialize the cache space
    call initialize_cache_space(negf,'G_r')
    call initialize_cache_space(negf,'G_n')

  end subroutine set_elph_block_dephasing

  !> Set values for the semi-local electron phonon dephasing model
  !! (elastic scattering only)
  subroutine set_elph_s_dephasing(negf, coupling, orbsperatom, niter)
    type(Tnegf), intent(inout) :: negf
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    integer,  dimension(:), allocatable, intent(in) :: orbsperatom
    integer, intent(in) :: niter

    type(TInteractionNode), pointer :: node
    call negf%interactList%add(node)
    call elphondephd_create(node%inter)
    select type(pInter => node%inter)
    type is(ElPhonDephS)
      call elphondephs_init(pInter, negf%str, coupling, orbsperatom, negf%S, niter)
    class default
      stop 'ERROR: error of type downcast to ElPhonDephS'
    end select

    ! Initialize the cache space
    call initialize_cache_space(negf,'G_r')
    call initialize_cache_space(negf,'G_n')

  end subroutine set_elph_s_dephasing

  subroutine set_elph_inelastic(negf, coupling, wq, Temp, dz, eps0, eps_inf, q0, area, niter)
    type(Tnegf), intent(inout) :: negf
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    real(dp), intent(in) :: wq
    real(dp), intent(in) :: Temp
    real(dp), intent(in) :: dz
    real(dp), intent(in) :: eps0 
    real(dp), intent(in) :: eps_inf
    real(dp), intent(in) :: q0
    real(dp), intent(in) :: area
    integer, intent(in) :: niter

    type(TInteractionNode), pointer :: node
    call negf%interactList%add(node)
    call elphondephd_create(node%inter)
    select type(pInter => node%inter)
    type is(ElPhonInel)
      call elphoninel_init(pInter, negf%cartComm%id, negf%str, negf%basis, coupling, &
          &  wq, Temp, dz, eps0, eps_inf, q0, area, niter)
    class default
      stop 'ERROR: error of type downcast to ElPhonInel'
    end select

    ! Initialize the cache space
    call initialize_cache_space(negf,'G_r')
    call initialize_cache_space(negf,'G_n')
    
  end subroutine set_elph_inelastic

  !> clean interactions objects
  subroutine destroy_interactions(this)
    class(Tnegf) :: this

    call this%interactList%destroy()

  end subroutine destroy_interactions

  subroutine initialize_cache_space(negf, matrix)
    type(Tnegf), intent(inout) :: negf 
    character(3), intent(in) :: matrix

    if (matrix == 'G_r') then
      if (.not.associated(negf%G_r)) then
         allocate(TMatrixCacheMem::negf%G_r)
         select type(p => negf%G_r)
         type is(TMatrixCacheMem)
            p%tagname='G_r'
         end select
      end if
    end if

    if (matrix == 'G_n') then
      if (.not.associated(negf%G_n)) then
         allocate(TMatrixCacheMem::negf%G_n)
         select type(p => negf%G_n)
         type is(TMatrixCacheMem)
            p%tagname='G_n'
         end select
      end if
    end if

  end subroutine initialize_cache_space

  subroutine set_defaults(this)
    class(Tnegf) :: this

     this%verbose = 10

     this%scratch_path = './GS/'
     this%out_path = './'
     this%DorE = 'D'           ! Density or En.Density

     this%ReadOldDM_SGFs = 1   ! Compute Surface G.F. do not save
     this%ReadOldT_SGFs = 1    ! Compute Surface G.F. do not save

     this%kwght = 1.0_dp
     this%ikpoint = 1

     this%Ec = 0.0_dp
     this%Ev = 0.0_dp
     this%DeltaEc = 0.0_dp
     this%DeltaEv = 0.0_dp

     this%E = 0.0_dp            ! Holding variable
     this%dos = 0.0_dp          ! Holding variable
     this%eneconv = 1.0_dp      ! Energy conversion factor

     this%isSid = .false.
     this%intHS = .true.
     this%intDM = .true.

     this%delta = 1.d-4      ! delta for G.F.
     this%dos_delta = 1.d-4  ! delta for DOS
     this%deltaModel = 1     ! deltaOmega model
     this%wmax = 0.009_dp     ! about 2000 cm^-1 cutoff
     this%Emin = 0.0_dp        ! Tunneling or dos interval
     this%Emax = 0.0_dp        !
     this%Estep = 0.0_dp       ! Tunneling or dos E step
     this%g_spin = 2.0_dp      ! spin degeneracy

     this%Np_n = (/20, 20/)  ! Number of points for n
     this%Np_p = (/20, 20/)  ! Number of points for p
     this%n_kt = 10          ! Numero di kT per l'integrazione
     this%n_poles = 3        ! Numero di poli
     this%activecont = 0     ! contact selfenergy
     this%min_or_max = 1     ! Set reference cont to max(mu)
     this%refcont = 1        ! call set_ref_cont()
     this%outer = 2          ! Compute full D.M. L,U extra
     this%dumpHS = .false.
     this%int_acc = 0.001_dp    ! Integration accuracy
                             ! Only in adaptive refinement
     this%ndos_proj = 0

     this%surface_green_cache = TMatrixCacheDisk(scratch_path=this%scratch_path)
     
     allocate(this%ni(1))
     allocate(this%nf(1))
     this%ni(1) = 1
     this%nf(1) = 2          !

   end subroutine set_defaults


   subroutine set_phph(negf, order, filename)
      type(Tnegf) :: negf
      integer, intent(in) :: order
      character(*), intent(in) :: filename

      print*,'(set_phph) init_phph'
      call init_phph(negf%phph, negf%str%central_dim, order, negf%str%mat_PL_start, negf%str%mat_PL_end)

      print*,'(set_phph) load_phph_coupl'
      call load_phph_couplings(negf%phph, filename)

   end subroutine set_phph



   subroutine print_all_vars(this,io)
     class(Tnegf) :: this
     integer, intent(in) :: io

     integer :: ii

     call print_Tstruct(this%str,io)

     write(io,*) 'verbose=',this%verbose

     write(io,*) 'scratch= "'//trim(this%scratch_path)//'"'
     write(io,*) 'output= "'//trim(this%out_path)//'"'

     write(io,*) 'isSid= ',this%isSid
     write(io,*) 'mu= ', this%mu
     write(io,*) 'mu_n= ', this%mu_n
     write(io,*) 'mu_p= ', this%mu_p
     write(io,*) 'T= ',this%kbT

     do ii = 1, this%str%num_conts
       write(io,*) 'Contact Parameters for ',trim(this%cont(ii)%name)
       write(io,*) 'mu= ', this%cont(ii)%mu
       write(io,*) 'WideBand= ', this%cont(ii)%FictCont
       write(io,*) 'DOS= ', this%cont(ii)%contact_DOS
       write(io,*) 'mu_n= ', this%cont(ii)%mu_n
       write(io,*) 'mu_p= ', this%cont(ii)%mu_p
       write(io,*) 'Density Matrix T= ',this%cont(ii)%kbT_dm
       write(io,*) 'Transmission T= ',this%cont(ii)%kbT_t
     end do

     write(io,*) 'Contour Parameters:'
     write(io,*) 'Ec= ', this%Ec
     write(io,*) 'Ev= ', this%Ev
     write(io,*) 'DEc= ', this%DeltaEc
     write(io,*) 'DEv= ', this%DeltaEv
     write(io,*) 'ReadoldDM_SGFs= ',this%ReadoldDM_SGFs
     write(io,*) 'ReadoldT_SGF= ',this%ReadoldT_SGFs
     write(io,*) 'Np_n= ',this%Np_n
     write(io,*) 'Np_p= ', this%Np_p
     write(io,*) 'nkT= ', this%n_kt
     write(io,*) 'nPoles= ', this%n_poles
     write(io,*) 'min_or_max= ', this%min_or_max
     write(io,*) 'refcont= ', this%refcont

     write(io,*) 'Transmission Parameters:'
     write(io,*) 'Emin= ',this%Emin
     write(io,*) 'Emax= ',this%Emax
     write(io,*) 'Estep= ',this%Estep
     write(io,*) 'g_spin= ',this%g_spin
     write(io,*) 'delta= ',this%delta
     write(io,*) 'dos_delta= ',this%dos_delta
     write(io,*) 'ni= ', this%ni
     write(io,*) 'nf= ', this%nf
     write(io,*) 'DorE= ',this%DorE

     write(io,*) 'Internal variables:'
     write(io,*) 'intHS= ',this%intHS
     write(io,*) 'intDM= ',this%intDM
     write(io,*) 'kp= ', this%ikpoint
     write(io,*) 'kwght= ', this%kwght
     write(io,*) 'E= ',this%E
     write(io,*) 'outer= ', this%outer
     write(io,*) 'DOS= ',this%dos
     write(io,*) 'Eneconv= ',this%eneconv
     write(io,*) 'activecont= ', this%activecont
  end subroutine print_all_vars

end module lib_param
